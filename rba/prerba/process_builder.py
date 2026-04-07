"""Module defining ProcessBuilder class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# local imports
import rba.xml
from rba.prerba.user_machinery import UserMachinery
from rba.prerba.macromolecule import Protein

class ProcessBuilder(object):
    """
    Class initializing defined process structure used by RBA.

    Attributes
    ----------
    data : rba.prerba.user_data.UserData
        User data.

    """

    def __init__(self,data):
        """
        Build object from user data.

        Parameters
        ----------
        data : rba.prerba.user_data.UserData
            User data.

        """
        self.data = data

    def create_processes_and_processing_maps(self):
        """Generates and returns Process and Processing map xml objects, based on process definition"""
        process_configuration=_extract_process_inputs_machinery_and_information(data=self.data)
        return(self.create_processes(process_configuration) , self.create_processing_maps(process_configuration) , process_configuration["Machinery"])


    def create_processing_maps(self,process_info):
        """Generate processing maps for processes from default processing maps"""
        processing_maps={}
        for process_id in process_info["Processes"].keys():
            template_processing_map_id=process_info["Processes"][process_id].get("TEMPLATE PROCESSING MAP","dummy_processing_map")
            processing_location=self.data.generate_mapping_rba_to_sbml_compartments(rba_compartment=process_info["Processes"][process_id].get("PROCESSING LOCATION","Cytoplasm"))
            macromolecule_type=process_info["Processes"][process_id]["MACROMOLECULE TYPE"]

            if process_info["Processes"][process_id]["Machinery"].has_nonempty_composition():
                has_machinery=True
            else:
                has_machinery=False

            processing_map_id=_get_processing_map_id(process_info,process_id,has_machinery)

            if template_processing_map_id == "default_processing":
                map=self.data.default.processing_maps.default_map(map_id=processing_map_id,macromolecule_type=macromolecule_type,non_spontaneous=has_machinery)
            elif template_processing_map_id == "default_processing_atp_dependent":
                map=self.data.default.processing_maps.default_map_atp_dependent(map_id=processing_map_id,processing_location=processing_location,macromolecule_type=macromolecule_type,non_spontaneous=has_machinery)            
            elif template_processing_map_id =="default_translation":
                map=self.data.default.processing_maps.default_translation_map(map_id=processing_map_id,processing_location=processing_location, cofactors=self.data.cofactors(), is_prokaryote=self.data._parameters.get("SPECIES_CATEGORY",'PROKARYOTE')=='PROKARYOTE',non_spontaneous=has_machinery)            
            elif template_processing_map_id =="default_transcription":
                map=self.data.default.processing_maps.default_transcription_map(map_id=processing_map_id,processing_location=processing_location,non_spontaneous=has_machinery) 
            elif template_processing_map_id =="default_dna_replication":
                map=self.data.default.processing_maps.default_replication_map(map_id=processing_map_id,processing_location=processing_location,non_spontaneous=has_machinery) 
            elif template_processing_map_id =="default_protein_degradation":
                map=self.data.default.processing_maps.default_protein_degradation_map(map_id=processing_map_id,processing_location=processing_location,cofactors=self.data.cofactors(),non_spontaneous=has_machinery) 
            elif template_processing_map_id =="default_rna_degradation":
                map=self.data.default.processing_maps.default_rna_degradation_map(map_id=processing_map_id,processing_location=processing_location,non_spontaneous=has_machinery) 
            else:
                map=self.data.default.processing_maps.default_map(map_id=processing_map_id,macromolecule_type=macromolecule_type,non_spontaneous=has_machinery)

            processing_maps[processing_map_id]=map
        
        return(list(processing_maps.values()))
    

    def create_processes(self,process_info):
        """Generate processes with machinery and processings"""
        processes=[]
        self.processing_input_fractions={}
        for process_id in process_info["Processes"].keys():
            process = rba.xml.Process("P_{}".format(process_id), process_id)
            # machinery
            if process_info["Processes"][process_id]["Machinery"]:
                for id_, sto in process_info["Processes"][process_id]["Machinery"].composition().items():
                    machine = process.machinery.machinery_composition
                    machine.reactants.append(rba.xml.SpeciesReference(id_, sto))

                capacity_param=process_info["Processes"][process_id].get("CAPACITY PARAMETER",float("nan"))
                if type(capacity_param) != str:
                    process.machinery.capacity.value = "chaperone_efficiency_CM"
                else:
                    process.machinery.capacity.value = capacity_param # get returns default chaperone_efficiency_LM if key is not present
            else:
                print("")
                print("WARNING: Process {} does not have a machinery composition, defined in fasta file ({}). \n " \
                "It is therefore defined as a machinery-less process.".format(process_id,self.data.input_path('machinery_composition/{}.fasta'.format(process_id))))

            process_has_no_inputs=True
            for mm_type in process_info["Processes"][process_id]["Inputs"].keys():
                if process_info["Processes"][process_id]["Inputs"][mm_type]:
                    process_has_no_inputs=False
                    if process_info["Processes"][process_id]["Machinery"].has_nonempty_composition():
                        has_machinery=True
                    else:
                        has_machinery=False
                    processing_map_id=_get_processing_map_id(process_info,process_id,has_machinery)

                    if process_info["Processes"][process_id].get("INPUT FRACTION",None) is not None:
                        input_fraction_parameter="{}__input_fraction".format(process_id)
                        self.processing_input_fractions[input_fraction_parameter]=process_info["Processes"][process_id].get("INPUT FRACTION",None)
                        processing = rba.xml.Processing(map_=processing_map_id, set_=mm_type, input_fraction_=input_fraction_parameter) #MMtype protein rna or dna
                    else:
                        processing = rba.xml.Processing(map_=processing_map_id, set_=mm_type) #MMtype protein rna or dna
                    mm_ids_in_inputs=[i.species for i in processing.inputs if type(i.species) == str]
                    for input_species in process_info["Processes"][process_id]["Inputs"][mm_type]:
                        if input_species not in mm_ids_in_inputs:
                            if type(input_species) == str:
                                processing.inputs.append(rba.xml.SpeciesReference(input_species, 1.0))
                                mm_ids_in_inputs.append(input_species)

                    if process_info["Processes"][process_id].get("PROCESSING TYPE","production")=="production":
                        process.processings.productions.append(processing)
                    if process_info["Processes"][process_id].get("PROCESSING TYPE","production")=="degradation":
                        process.processings.degradations.append(processing)        
            if process_has_no_inputs:
                print("")
                print("WARNING: Process {} does not have any assigned inputs. (Consider removing it?)".format(process_id))
            processes.append(process)

        return processes        
        

def _get_processing_map_id(process_info,process_id,has_machinery=True):
    """Create ID of processing map, based on process information"""
    template_processing_map_id=process_info["Processes"][process_id].get("TEMPLATE PROCESSING MAP","dummy_processing_map")
    processing_location=process_info["Processes"][process_id].get("PROCESSING LOCATION","Cytoplasm")
    macromolecule_type=process_info["Processes"][process_id]["MACROMOLECULE TYPE"]
    if not has_machinery:
        suffix="__spontaneous"
    else:
        suffix=""

    if template_processing_map_id == "default_processing":
        processing_map_id="{}__{}__{}{}".format(template_processing_map_id,
                                            macromolecule_type,
                                            processing_location,
                                            suffix)
    elif template_processing_map_id == "default_processing_atp_dependent":
        processing_map_id="{}__{}__{}{}".format(template_processing_map_id,
                                            macromolecule_type,
                                            processing_location,
                                            suffix)
                
    elif template_processing_map_id =="default_translation":
        processing_map_id="{}__{}{}".format("translation",
                                        processing_location,
                                        suffix)
                
    elif template_processing_map_id =="default_transcription":
        processing_map_id="{}__{}{}".format("transcription",
                                        processing_location,
                                        suffix)
                
    elif template_processing_map_id =="default_dna_replication":
        processing_map_id="{}__{}{}".format("dna_replication",
                                        processing_location,
                                        suffix)
                
    elif template_processing_map_id =="default_protein_degradation":
        processing_map_id="{}__{}{}".format("protein_degradation",
                                        processing_location,
                                        suffix)
                
    elif template_processing_map_id =="default_rna_degradation":
        processing_map_id="{}__{}{}".format("rna_degradation",
                                        processing_location,
                                        suffix)
                
    else:
        processing_map_id="{}__{}__{}{}".format("default_processing",
                                            macromolecule_type,
                                            processing_location,
                                            suffix)
    return(processing_map_id)


def _extract_process_inputs_machinery_and_information(data):
    """Construct sets of processing inputs and machinery components for process from process definition and machinery composition info"""
    # lists of currently defined macromolecules as potential inputs to processes 
    # (process machinery subunit macromolecules from fasta files are not defined yet)
    all_inputs={"protein":list(data.enzymatic_localised_proteins+_build_dummy_average_proteins(data)), #enzymatic and average proteins
                      "rna":list(list(data.trnas)+list(data.default_mrnas.values())), #trnas and default mrnas
                      "dna": list(data.default_dnas.values()) # default dnas
                      }
    
    # inputs for each process, resolved by macromolecule type
    inputs={row["PROCESS"]:{"protein":[],"rna":[],"dna":[]} for row in data.process_definition}

    machinery={}
    process_info={}

    #assign inputs to processes and add subunits of machinery as inputs iteratively 
    # until no new processes are required anymore
    macromolecules_to_be_assigned=True
    while macromolecules_to_be_assigned:
        added_processes=[]
        for defined_process in data.process_definition: # iterate over entries in process definition (rows in file)
            
            process_id=defined_process["PROCESS"]

            process_info[process_id]={}
            process_info[process_id].update(defined_process)

            macromolecule_specific_id=defined_process.get("ID",None)
            macromolecule_location=defined_process.get("LOCATION",None)
            macromolecule_origin=defined_process.get("ORIGIN",None)
            macromolecule_type=defined_process.get("MACROMOLECULE TYPE",None)
            
            ### Define inputs to process, based on macromolecule type, id, location and origin
            if type(macromolecule_specific_id) == str:
                #add macromolecule, specified by ID, as input to this process (if not already in list)
                inputs[process_id][macromolecule_type] += [i for i in [macromolecule_specific_id] if i not in inputs[process_id][macromolecule_type]]
            else:
                if type(macromolecule_origin) == str:
                    if type(macromolecule_location) == str:
                        #add macromolecules, with respective ORIGIN AND LOCATION, as input to this process (if not already in list)
                        inputs[process_id][macromolecule_type] += [mm.id for mm in all_inputs[macromolecule_type] if mm.origin==macromolecule_origin and mm.location==macromolecule_location and mm.id not in inputs[process_id][macromolecule_type]]
                    else:
                        #add macromolecules, with respective ORIGIN, as input to this process (if not already in list)
                        inputs[process_id][macromolecule_type] += [mm.id for mm in all_inputs[macromolecule_type] if mm.origin==macromolecule_origin and mm.id not in inputs[process_id][macromolecule_type]]
                else:
                    if type(macromolecule_location) == str:
                        #add macromolecules, with respective LOCATION, as input to this process (if not already in list)
                        inputs[process_id][macromolecule_type] += [mm.id for mm in all_inputs[macromolecule_type] if mm.location==macromolecule_location and mm.id not in inputs[process_id][macromolecule_type]]
                    else:
                        #add macromolecules of specific TYPE, as input to this process (if not already in list)
                        inputs[process_id][macromolecule_type] += [mm.id for mm in all_inputs[macromolecule_type] if mm.id not in inputs[process_id][macromolecule_type]]

            ### define machinery of processes not already accounted for

            # check if process-machinery has not already been defined previously
            if process_id not in machinery:
                added_processes.append(process_id) # mark process as added in this iteration
                # define process machinery from fata files
                machinery[process_id] = UserMachinery(filename=data.input_path('machinery_composition/{}.fasta'.format(process_id)),
                                                      protein_data=data.protein_data,
                                                      location_separator=data._macromolecule_location_separator,
                                                      metabolic_species=[sbml_species.id for sbml_species in data.sbml_data.species])

        ### If there were processes added in this iteration
        if added_processes:
            # Define their subunits as input macromolecules for next iteration
            all_inputs={"protein":[],"rna":[],"dna":[]}
            for process in added_processes:
                all_inputs["protein"]+=[i for i in machinery[process].proteins if i not in all_inputs["protein"]]
                all_inputs["rna"]+=[i for i in machinery[process].rnas if i not in all_inputs["rna"]]
        # If no processes and therefore machinery subunits were added in this iteration
        else:
            macromolecules_to_be_assigned=False
            # --> STOPPING criterion

    ### Generate output
    out={"Processes":{},"Machinery":{"Proteins":[],"RNAs":[]}}
    for process in list(set(list(list(machinery.keys())+list(inputs.keys())))):
        out["Processes"][process]={}
        out["Processes"][process].update(process_info[process])
        out["Processes"][process]["Inputs"]=inputs[process]

        if process in machinery: # does process have machinery
            # if so, define it
            out["Processes"][process]["Machinery"]=machinery[process]
            out["Machinery"]["Proteins"]+=machinery[process].proteins
            out["Machinery"]["RNAs"]+=machinery[process].rnas
        else:
            out["Processes"][process]["Machinery"]=None

    return(out)


def _build_dummy_average_proteins(data):
    """Build dummy of average protein as proxy for process"""
    out=[]
    for compartment in data.compartments():
        protein = Protein()
        protein.id = data.protein_data.average_protein_id(compartment)
        protein.location = compartment
        protein.origin = data.default_genome_location
        out.append(protein)
    return(out)


