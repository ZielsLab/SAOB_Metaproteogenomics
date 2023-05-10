###SECTION 1: Import dependencies and create blank model

from cobra import Model, Reaction, Metabolite
import cobra
import xlsxwriter

model = Model('model')

#Extracellular metabolites
h2o_e = Metabolite('h2o_e', formula='H2O', name='H2O', compartment='e', charge=0)
pi_e = Metabolite('pi_e', formula='HO4P', name='Phosphate', compartment='e', charge=-2)
co2_e = Metabolite('co2_e', formula='CO2', name='CO2', compartment='e', charge=0)
for_e = Metabolite('for_e', formula='CHO2', name='Formate', compartment='e', charge=-1)
h_e = Metabolite('h_e', formula='H', name='H+', compartment='e', charge=1)
ac_e = Metabolite('ac_e', formula='C2H3O2', name='Acetate', compartment='e', charge=-1)
h2_e = Metabolite('h2_e', formula='H2', name='Hydrogen', compartment='e', charge=0)
ch4_e = Metabolite('ch4_e', formula='CH4', name='Methane', compartment='e', charge=0)

#Counting metabolite
atp_COMM =Metabolite('atp_COMM', formula='', name='', compartment='e', charge=0)

#DTU068_1 metabolites
h2o_DTU068_1e = Metabolite('h2o_DTU068_1e', formula='H2O', name='H2O', compartment='e', charge=0)
co2_DTU068_1e = Metabolite('co2_DTU068_1e', formula='CO2', name='CO2', compartment='e', charge=0)
for_DTU068_1e = Metabolite('for_DTU068_1e', formula='CHO2', name='Formate', compartment='e', charge=-1)
h_DTU068_1e = Metabolite('h_DTU068_1e', formula='H', name='H+', compartment='e', charge=1)
ac_DTU068_1e = Metabolite('ac_DTU068_1e', formula='C2H3O2', name='Acetate', compartment='e', charge=-1)
h2_DTU068_1e = Metabolite('h2_DTU068_1e', formula='H2', name='Hydrogen', compartment='e', charge=0)

ac_DTU068_1c = Metabolite('ac_DTU068_1c', formula='C2H3O2', name='Acetate', compartment='DTU068_1c', charge=-1)
atp_DTU068_1c = Metabolite('atp_DTU068_1c', formula='C10H12N5O13P3', name='ATP', compartment='DTU068_1c', charge=-4)
adp_DTU068_1c = Metabolite('adp_DTU068_1c', formula='C10H12N5O10P2', name='ADP', compartment='DTU068_1c', charge=-3)
h_DTU068_1c = Metabolite('h_DTU068_1c', formula='H', name='H+', compartment='DTU068_1c', charge=1)
pi_DTU068_1c = Metabolite('pi_DTU068_1c', formula='HO4P', name='xylose-D', compartment='DTU068_1c', charge=-2)
h2o_DTU068_1c = Metabolite('h2o_DTU068_1c', formula='H2O', name='H2O', compartment='DTU068_1c', charge=0)
actp_DTU068_1c = Metabolite('actp_DTU068_1c', formula='C2H3O5P', name='Acetyl phosphate', compartment='DTU068_1c', charge=-2)
for_DTU068_1c = Metabolite('for_DTU068_1c', formula='CHO2', name='Formate', compartment='DTU068_1c', charge= -1)
accoa_DTU068_1c = Metabolite('accoa_DTU068_1c', formula='C23H34N7O17P3S', name='Acetyl-CoA', compartment='DTU068_1c', charge=-4)
coa_DTU068_1c = Metabolite('coa_DTU068_1c', formula='C21H32N7O16P3S', name='Coenzyme A', compartment='DTU068_1c', charge=-4)
co_DTU068_1c = Metabolite('co_DTU068_1c', formula='CO', name='Carbon monoxide', compartment='DTU068_1c', charge=0)
fdred_DTU068_1c = Metabolite('fdred_DTU068_1c', formula='Fe8S8X', name='Ferredoxin (reduced) 2[4Fe-4S]', compartment='DTU068_1c', charge= -2)
fdox_DTU068_1c = Metabolite('fdox_DTU068_1c', formula='Fe8S8X', name='Ferredoxin (oxidized) 2[4Fe-4S]', compartment='DTU068_1c', charge= 0)
co2_DTU068_1c = Metabolite('co2_DTU068_1c', formula='CO2', name='CO2', compartment='DTU068_1c', charge= 0)
cfesp_DTU068_1c = Metabolite('cfesp_DTU068_1c', formula='C19CoN4R21', name='Corrinoid Iron sulfur protein', compartment='DTU068_1c', charge=-1)
mecfsp_DTU068_1c = Metabolite('mecfsp_DTU068_1c', formula='C20H3CoN4R21', name='Methylcorrinoid iron sulfur protein', compartment='DTU068_1c', charge=0)
_5mthf_DTU068_1c = Metabolite('_5mthf_DTU068_1c', formula='C20H24N7O6', name='5-Methyltetrahydrofolate', compartment='DTU068_1c', charge=-1)
mlthf_DTU068_1c = Metabolite('mlthf_DTU068_1c', formula='C20H21N7O6', name='5,10-Methylenetetrahydrofolate', compartment='DTU068_1c', charge=-2)
methf_DTU068_1c = Metabolite('methf_DTU068_1c', formula='C20H20N7O6', name='5,10-Methenyltetrahydrofolate', compartment='DTU068_1c', charge=-1)
thf_DTU068_1c = Metabolite('thf_DTU068_1c', formula='C19H21N7O6', name='5,6,7,8-Tetrahydrofolate', compartment='DTU068_1c', charge=-2)
_10fthf_DTU068_1c = Metabolite('_10fthf_DTU068_1c', formula='C20H21N7O7', name='10-Formyltetrahydrofolate', compartment='DTU068_1c', charge=-2)
h2_DTU068_1c = Metabolite('h2_DTU068_1c', formula='H2', name='Hydrogen', compartment='DTU068_1c', charge= 0)
nad_DTU068_1c = Metabolite('nad_DTU068_1c', formula='C21H26N7O14P2', name='Nicotinamide adenine dinucleotide', compartment='DTU068_1c', charge=-1)
nadh_DTU068_1c = Metabolite('nadh_DTU068_1c', formula='C21H27N7O14P2', name='Nicotinamide adenine dinucleotide - reduced', compartment='DTU068_1c', charge=-2)
h_DTU068_1i = Metabolite('h_DTU068_1i', formula='H', name='H+', compartment='DTU068_1i', charge=1)


reaction = Reaction('DTU068_1_ACK')
reaction.name = 'DTU068_1: Acetate kinase'
reaction.subsystem = 'Acetate Metabolism'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({ac_DTU068_1c: -1.0,
                         atp_DTU068_1c: -1.0,
                         actp_DTU068_1c: 1.0,
                         adp_DTU068_1c: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###PTA###

reaction = Reaction('DTU068_1_PTA')
reaction.name = 'DTU068_1: Phosphotransacetylase'
reaction.subsystem = 'Acetate Metabolism'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({accoa_DTU068_1c: 1.0,
                          pi_DTU068_1c: 1.0,
                          actp_DTU068_1c: -1.0,
                          coa_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###CODH/ACS###

#co_DTU068_1c + coa_DTU068_1c  + mecfsp_DTU068_1c -> accoa_DTU068_1c + cfesp_DTU068_1c + h_DTU068_1c

reaction = Reaction('DTU068_1_ACS')
reaction.name = 'DTU068_1: Acetyl-CoA synthase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co_DTU068_1c: 1.0,
                            coa_DTU068_1c: 1.0,
                            mecfsp_DTU068_1c: 1.0,
                            accoa_DTU068_1c: -1.0,
                            cfesp_DTU068_1c: -1.0,
                            h_DTU068_1c: -1.0})

model.add_reactions([reaction])

#co2_DTU068_1c + 2.0 h_DTU068_1c + fdred_DTU068_1c <-> h2o_DTU068_1c + co_DTU068_1c + fdox__DTU068_1c

#BIGG uses a different form of ferredoxin

reaction = Reaction('DTU068_1_CODH')
reaction.name = 'DTU068_1: Carbon monoxide dehydrogenase / acetyl-CoA synthase 2'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({co2_DTU068_1c: 1.0,
                            h_DTU068_1c: 2.0,
                            fdred_DTU068_1c: 1.0,
                            h2o_DTU068_1c: -1.0,
                            co_DTU068_1c: -1.0,
                            fdox_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


#5mthf_DTU068_1c + cfesp_DTU068_1c -> thf_DTU068_1c + mecfsp_DTU068_1c

cfesp_DTU068_1c = Metabolite('cfesp_DTU068_1c', formula='C19CoN4R21', name='Corrinoid Iron sulfur protein', compartment='DTU068_1c', charge=-1)
mecfsp_DTU068_1c = Metabolite('mecfsp_DTU068_1c', formula='C20H3CoN4R21', name='Methylcorrinoid iron sulfur protein', compartment='DTU068_1c', charge=0)

reaction = Reaction('DTU068_1_AcsE')
reaction.name = 'DTU068_1: Methyltetrahydrofolate:corrinoid/iron-sulfur protein methyltransferase (MeTr)'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_5mthf_DTU068_1c: 1.0,
                            cfesp_DTU068_1c: 1.0,
                            thf_DTU068_1c: -1.0,
                            mecfsp_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###MeTF (MTHFR)###
#Electron flow from Methyl-THF to MvhD to HdrABC to MQ to cyt-b to FdxH to FdhA
#Consumes two protons from ion motive force

reaction = Reaction('DTU068_1_MetF_MetV')
reaction.name = 'DTU068_1: MetF + MetV'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_5mthf_DTU068_1c: -1.0,
                            co2_e: -1.0,
                            h_DTU068_1i: -2.0,
                            h_e: -1.0,
                            mlthf_DTU068_1c: 1.0,
                            h_DTU068_1c: 5.0,
                            for_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

###MeTF (MTHFR)###
#Electron flow from Methyl-THF to MvhD to HdrABC to MQ to cyt-b to HysAB
#Consumes two protons from ion motive force

reaction = Reaction('DTU068_1_HysAB')
reaction.name = 'DTU068_1: HysAB'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_5mthf_DTU068_1c: -1.0,
                            h_DTU068_1i: -2.0,
                            h_e: -2.0,
                            mlthf_DTU068_1c: 1.0,
                            h_DTU068_1c: 5.0,
                            h2_e: 1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_FoID_MTHFD_1')
reaction.name = 'DTU068_1: Methylenetetrahydrofolate dehydrogenase NAD'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({methf_DTU068_1c: 1.0,
                            nadh_DTU068_1c: 1.0,
                            mlthf_DTU068_1c: -1.0,
                            nad_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_FoID_MTHFD_2')
reaction.name = 'DTU068_1: Methenyltetrahydrofolate cyclohydrolase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({_10fthf_DTU068_1c: 1.0,
                            h_DTU068_1c: 1.0,
                            h2o_DTU068_1c: -1.0,
                            methf_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_Fhs(FTHFL)')
reaction.name = 'DTU068_1: Formate-tetrahydrofolate ligase'
reaction.subsystem = 'Wood Ljungadhl Pathway'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({for_DTU068_1c: 1.0,
                            atp_DTU068_1c: 1.0,
                            thf_DTU068_1c: 1.0,
                            _10fthf_DTU068_1c: -1.0,
                            adp_DTU068_1c: -1.0,
                            pi_DTU068_1c: -1.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_Fdh_NuoEF')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'DTU068_1: Electron Confurcating Formate Dehydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({nadh_DTU068_1c: 1.0,
                          fdred_DTU068_1c: 1.0,
                          h_DTU068_1c: 1.0,
                          co2_DTU068_1c: 2.0,
                          nad_DTU068_1c: -1.0,
                          fdox_DTU068_1c: -1.0,
                          for_DTU068_1c: -2.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))


###Hydrogen Production###

reaction = Reaction('DTU068_1_ECH')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'DTU068_1: Energy conserving hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_DTU068_1c: -1.0,
                          h_DTU068_1c: -4.0,
                          h2_DTU068_1c: 1.0,
                          fdox_DTU068_1c: 1.0,
                          h_DTU068_1i: 2.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_HYDABC')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'DTU068_1: HydABC'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_DTU068_1c: -1.0,
                          nadh_DTU068_1c: -1.0,
                          h_DTU068_1c: -3.0,
                          h2_DTU068_1c: 2.0,
                          nad_DTU068_1c: 1.0,
                          fdox_DTU068_1c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_NiFeH2ASE')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'DTU068_1: NiFe Hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({nadh_DTU068_1c: -1.0,
                          h_DTU068_1c: -1.0,
                          h2_DTU068_1c: 1.0,
                          nad_DTU068_1c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_ATPS4r')
#This reaction differs from the BiGG reaction because this model assumes a different compartment for ion motive force generation
reaction.name = 'DTU068_1: ATP Synthase'
reaction.subsystem = 'Energy Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({adp_DTU068_1c: -1.0,
                          pi_DTU068_1c: -1.0,
                          h_DTU068_1i: -4.0,
                          atp_DTU068_1c: 1.0,
                          h_DTU068_1c: 3.0,
                          h2o_DTU068_1c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ATP Hydrolysis

#atp_DTU068_1c + h2o_DTU068_1c <-> adp_DTU068_1c + pi_DTU068_1c + h_DTU068_1c + ATP_COMM_e

reaction = Reaction('DTU068_1_ATP_Hydrolysis')
reaction.name = 'DTU068_1: ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_DTU068_1c: -1.0,
                          h2o_DTU068_1c: -1.0,
                          adp_DTU068_1c: 1.0,
                          pi_DTU068_1c: 1.0,
                          h_DTU068_1c: 1.0,
                          atp_COMM: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_EX_H2')
reaction.name = 'DTU068_1: EX - H2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2_DTU068_1c: 1.0,
                          h2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_EX_Acetate')
reaction.name = 'DTU068_1: Transport - Acetate'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ac_DTU068_1c: 1.0,
                          h_DTU068_1c: 1.0,
                          ac_e: -1.0,
                          h_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_EX_Protons')
reaction.name = 'DTU068_1: Transport - Protons'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h_DTU068_1c: 1.0,
                          h_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('DTU068_1_EX_h2o')
reaction.name = 'DTU068_1: Transport - h2o'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2o_DTU068_1c: 1.0,
                          h2o_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('DTU068_1_EX_co2')
reaction.name = 'DTU068_1: Transport - co2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({co2_DTU068_1c: 1.0,
                          co2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('DTU068_1_EX_for')
reaction.name = 'DTU068_1: Transport - for'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({for_DTU068_1c: 1.0,
                          h_DTU068_1c: 1.0,
                          h_e: -1.0,
                          for_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


#Met_1 metabolites
h2o_Met_1e = Metabolite('h2o_Met_1e', formula='H2O', name='H2O', compartment='e', charge=0)
co2_Met_1e = Metabolite('co2_Met_1e', formula='CO2', name='CO2', compartment='e', charge=0)
for_Met_1e = Metabolite('for_Met_1e', formula='CHO2', name='Formate', compartment='e', charge=-1)
h_Met_1e = Metabolite('h_Met_1e', formula='H', name='H+', compartment='e', charge=1)
h2_Met_1e = Metabolite('h2_Met_1e', formula='H2', name='Hydrogen', compartment='e', charge=0)
ch4_Met_1e = Metabolite('ch4_Met_1e', formula='CH4', name='Methane', compartment='e', charge=0)

ac_Met_1c = Metabolite('ac_Met_1c', formula='C2H3O2', name='Acetate', compartment='Met_1c', charge=-1)
atp_Met_1c = Metabolite('atp_Met_1c', formula='C10H12N5O13P3', name='ATP', compartment='Met_1c', charge=-4)
adp_Met_1c = Metabolite('adp_Met_1c', formula='C10H12N5O10P2', name='ADP', compartment='Met_1c', charge=-3)
h_Met_1c = Metabolite('h_Met_1c', formula='H', name='H+', compartment='Met_1c', charge=1)
pi_Met_1c = Metabolite('pi_Met_1c', formula='HO4P', name='xylose-D', compartment='Met_1c', charge=-2)
h2o_Met_1c = Metabolite('h2o_Met_1c', formula='H2O', name='H2O', compartment='Met_1c', charge=0)
actp_Met_1c = Metabolite('actp_Met_1c', formula='C2H3O5P', name='Acetyl phosphate', compartment='Met_1c', charge=-2)
for_Met_1c = Metabolite('for_Met_1c', formula='CHO2', name='Formate', compartment='Met_1c', charge= -1)
accoa_Met_1c = Metabolite('accoa_Met_1c', formula='C23H34N7O17P3S', name='Acetyl-CoA', compartment='Met_1c', charge=-4)
coa_Met_1c = Metabolite('coa_Met_1c', formula='C21H32N7O16P3S', name='Coenzyme A', compartment='Met_1c', charge=-4)
co_Met_1c = Metabolite('co_Met_1c', formula='CO', name='Carbon monoxide', compartment='Met_1c', charge=0)
fdred_Met_1c = Metabolite('fdred_Met_1c', formula='Fe8S8X', name='Ferredoxin (reduced) 2[4Fe-4S]', compartment='Met_1c', charge= -2)
fdox_Met_1c = Metabolite('fdox_Met_1c', formula='Fe8S8X', name='Ferredoxin (oxidized) 2[4Fe-4S]', compartment='Met_1c', charge= 0)
co2_Met_1c = Metabolite('co2_Met_1c', formula='CO2', name='CO2', compartment='Met_1c', charge= 0)
cfesp_Met_1c = Metabolite('cfesp_Met_1c', formula='C19CoN4R21', name='Corrinoid Iron sulfur protein', compartment='Met_1c', charge=-1)
mecfsp_Met_1c = Metabolite('mecfsp_Met_1c', formula='C20H3CoN4R21', name='Methylcorrinoid iron sulfur protein', compartment='Met_1c', charge=0)
_5mthf_Met_1c = Metabolite('_5mthf_Met_1c', formula='C20H24N7O6', name='5-Methyltetrahydrofolate', compartment='Met_1c', charge=-1)
mlthf_Met_1c = Metabolite('mlthf_Met_1c', formula='C20H21N7O6', name='5,10-Methylenetetrahydrofolate', compartment='Met_1c', charge=-2)
methf_Met_1c = Metabolite('methf_Met_1c', formula='C20H20N7O6', name='5,10-Methenyltetrahydrofolate', compartment='Met_1c', charge=-1)
thf_Met_1c = Metabolite('thf_Met_1c', formula='C19H21N7O6', name='5,6,7,8-Tetrahydrofolate', compartment='Met_1c', charge=-2)
_10fthf_Met_1c = Metabolite('_10fthf_Met_1c', formula='C20H21N7O7', name='10-Formyltetrahydrofolate', compartment='Met_1c', charge=-2)
h2_Met_1c = Metabolite('h2_Met_1c', formula='H2', name='Hydrogen', compartment='Met_1c', charge= 0)
nad_Met_1c = Metabolite('nad_Met_1c', formula='C21H26N7O14P2', name='Nicotinamide adenine dinucleotide', compartment='Met_1c', charge=-1)
nadh_Met_1c = Metabolite('nadh_Met_1c', formula='C21H27N7O14P2', name='Nicotinamide adenine dinucleotide - reduced', compartment='Met_1c', charge=-2)
h_Met_1i = Metabolite('h_Met_1i', formula='H', name='H+', compartment='Met_1i', charge=1)
fdxo_42_Met_1c = Metabolite('fdxo_42_Met_1c', formula='Fe12S12X', name='Ferredoxin - oxidized', compartment='Met_1c',charge=0)
fdxr_42_Met_1c = Metabolite('fdxr_42_Met_1c', formula='Fe12S12X', name='Ferredoxin - reduced', compartment='Met_1c',charge=-3)

mfr_b_Met_1c = Metabolite('mfr_b_Met_1c', formula='C34H44N6O14R', name='Methanofuran b', compartment='Met_1c', charge=-1)
formmfr_b_Met_1c = Metabolite('formmfr_b_Met_1c', formula='C35H43N6O15R', name='Formylmethanofuran b', compartment='Met_1c',charge=-2)
formh4spt_Met_1c = Metabolite('formh4spt_Met_1c', formula='C36H48N7O20P1', name='Formyltetrahydrosarcinapterin', compartment='Met_1c', charge=0)
h4spt_Met_1c = Metabolite('h4spt_Met_1c', formula='C35H48N7O19P1', name='Tetrahydrosarcinapterin', compartment='Met_1c',charge=0)
menylh4spt_Met_1c = Metabolite('menylh4spt_Met_1c', formula='C36H47N7O19P1', name='Methenyl-tetrahydrosarcinapterin',compartment='Met_1c', charge=1)
f420_2_Met_1c = Metabolite('f420_2_Met_1c', name='Coenzyme ferredoxin 420-2 (oxidized)', formula='C29H34N5O18P',compartment='Met_1c', charge=0)
f420_2h2_Met_1c = Metabolite('f420_2h2_Met_1c', name='Coenzyme ferredoxin 420-2 (reduced)', formula='C29H36O18N5P1',compartment='Met_1c', charge=0)
mleneh4spt_Met_1c = Metabolite('mleneh4spt_Met_1c', name='N5,N10-methylee-5,6,7,8-tetrahydromethanopterin', formula='C36H48N7O19P1', compartment='Met_1c', charge=0)
com_Met_1c = Metabolite('com_Met_1c', name='Coenzyme m', formula='C2H5O3S2', compartment='Met_1c', charge=0)
mcom_Met_1c = Metabolite('mcom_Met_1c', name='Methylcoenzyme m', formula='C3O3S2H7', compartment='Met_1c', charge=0)
mphen_Met_1c = Metabolite('mphen_Met_1c', name='Methanophenazine (oxidized)', formula='C37N2O1H50', compartment='Met_1c',charge=0)
mphenh2_Met_1c = Metabolite('mphenh2_Met_1c', name='Methanophenazine (reduced)', formula='C37N2O1H52', compartment='Met_1c',charge=0)
mh4spt_Met_1c = Metabolite('mh4spt_Met_1', name='N5-methyl-tetrahydrosarcinapterin', formula='C36H50N7O19P1',compartment='Met_1c', charge=0)

ch4_Met_1c = Metabolite('ch4_Met_1c', formula='CH4', name='Methane', compartment='Met_1c', charge=0)
cob_Met_1c = Metabolite('cob_Met_1c', formula='C11H19N1O7P1S1', name='Conenzyme B', compartment='Met_1c', charge=0)
hsfd_Met_1c = Metabolite('hsfd_Met_1c', formula='C13H22N1O10P1S3', name='Heterodisulfide', compartment='Met_1c', charge=0)


# [1] Formylmethanofuran Dehydrogenase - FMFD_b
reaction = Reaction('Met_1_FwdA-H')  # What do we do with R's if they're not balancing?
reaction.name = 'Met_1: Formylmethanofuran Dehydrogenase'
reaction.subsystem = 'Met_1'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({co2_Met_1c: -1.0,
                          fdred_Met_1c: -1.0,
                          h_Met_1c: -1.0,
                          mfr_b_Met_1c: -1.0,
                          fdox_Met_1c: 1.0,
                          h2o_Met_1c: 1.0,
                          formmfr_b_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# [2] Formylmethanofuran-tetrahydromethanopterin N-formyltransferase - FMFTSPFT_b
reaction = Reaction('Met_1_Ftr')
reaction.name = 'Met_1: Formylmethanofuran-Tetrahydromethanopterin N-Formyltransferase'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_1c: -1.0,
                          h4spt_Met_1c: -1.0,
                          mfr_b_Met_1c: 1.0,
                          formmfr_b_Met_1c: -1.0,
                          formh4spt_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# [3]? MTSPC
reaction = Reaction('Met_1_Mch')
reaction.name = 'Met_1: Mch'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_1c: -1.0,
                          h2o_Met_1c: 1.0,
                          menylh4spt_Met_1c: 1.0,
                          formh4spt_Met_1c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4RHi [9]
reaction = Reaction('Met_1_FrhABG')
reaction.name = 'Met_1: FrhABG'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({h2_Met_1c: -1.0,
                          f420_2_Met_1c: -1.0,
                          f420_2h2_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4MTSPD [4]
reaction = Reaction('Met_1_Mtd')
reaction.name = 'Met_1: Mtd'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_1c: 1.0,
                          f420_2_Met_1c: 1.0,
                          f420_2h2_Met_1c: -1.0,
                          mleneh4spt_Met_1c: 1.0,
                          menylh4spt_Met_1c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4MTSPR [5]
reaction = Reaction('Met_1_Mer')
reaction.name = 'Met_1: Mer'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({f420_2_Met_1c: 1.0,
                          mh4spt_Met_1c: 1.0,
                          f420_2h2_Met_1c: -1.0,
                          mleneh4spt_Met_1c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# MTSPCMMT [6]
reaction = Reaction('Met_1_MtrA-H')
reaction.name = 'Met_1: MtrA-H'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_1c: -2.0,
                          h_Met_1i: 2.0,
                          h4spt_Met_1c: 1.0,
                          mh4spt_Met_1c: -1.0,
                          com_Met_1c: -1.0,
                          mcom_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# MCR
reaction = Reaction('Met_1_McrABCDG')
reaction.name = 'Met_1: McrABCDG'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ch4_Met_1c: 1.0,
                          cob_Met_1c: -1.0,
                          mcom_Met_1c: -1.0,
                          hsfd_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# HDR
reaction = Reaction('Met_1_MvhADG+HdrABC')
reaction.name = 'Met_1: HDR'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({h2_Met_1c: -2.0,
                          h_Met_1c: 2.0,
                          hsfd_Met_1c: -1.0,
                          cob_Met_1c: 1.0,
                          com_Met_1c: 1.0,
                          fdox_Met_1c: -1.0,
                          fdred_Met_1c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

###Hydrogen Production###

reaction = Reaction('Met_1_ECH')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'Met_1: Energy conserving hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_Met_1c: -1.0,
                          h_Met_1c: -4.0,
                          h2_Met_1c: 1.0,
                          fdox_Met_1c: 1.0,
                          h_Met_1i: 2.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_1_ATPS4r')
#This reaction differs from the BiGG reaction because this model assumes a different compartment for ion motive force generation
reaction.name = 'Met_1: ATP Synthase'
reaction.subsystem = 'Energy Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({adp_Met_1c: -1.0,
                          pi_Met_1c: -1.0,
                          h_Met_1i: -4.0,
                          atp_Met_1c: 1.0,
                          h_Met_1c: 3.0,
                          h2o_Met_1c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ATP Hydrolysis

#atp_Met_1c + h2o_Met_1c <-> adp_Met_1c + pi_Met_1c + h_Met_1c + ATP_COMM_e

reaction = Reaction('Met_1_ATP_Hydrolysis')
reaction.name = 'Met_1: ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_Met_1c: -1.0,
                          h2o_Met_1c: -1.0,
                          adp_Met_1c: 1.0,
                          pi_Met_1c: 1.0,
                          h_Met_1c: 1.0,
                          atp_COMM: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('Met_1_EX_H2')
reaction.name = 'Met_1: EX - H2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2_Met_1c: 1.0,
                          h2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_1_EX_Protons')
reaction.name = 'Met_1: Transport - Protons'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h_Met_1c: 1.0,
                          h_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_1_EX_ch4')
reaction.name = 'Met_1: EX - ch4'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ch4_Met_1c: 1.0,
                          ch4_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_1_EX_h2o')
reaction.name = 'Met_1: EX - h2o'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2o_Met_1c: 1.0,
                          h2o_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_1_EX_co2')
reaction.name = 'Met_1: EX - co2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({co2_Met_1c: 1.0,
                          co2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#Met_2 metabolites
h2o_Met_2e = Metabolite('h2o_Met_2e', formula='H2O', name='H2O', compartment='e', charge=0)
co2_Met_2e = Metabolite('co2_Met_2e', formula='CO2', name='CO2', compartment='e', charge=0)
for_Met_2e = Metabolite('for_Met_2e', formula='CHO2', name='Formate', compartment='e', charge=-1)
h_Met_2e = Metabolite('h_Met_2e', formula='H', name='H+', compartment='e', charge=1)
h2_Met_2e = Metabolite('h2_Met_2e', formula='H2', name='Hydrogen', compartment='e', charge=0)
ch4_Met_2e = Metabolite('ch4_Met_2e', formula='CH4', name='Methane', compartment='e', charge=0)

ac_Met_2c = Metabolite('ac_Met_2c', formula='C2H3O2', name='Acetate', compartment='Met_2c', charge=-1)
atp_Met_2c = Metabolite('atp_Met_2c', formula='C10H12N5O13P3', name='ATP', compartment='Met_2c', charge=-4)
adp_Met_2c = Metabolite('adp_Met_2c', formula='C10H12N5O10P2', name='ADP', compartment='Met_2c', charge=-3)
h_Met_2c = Metabolite('h_Met_2c', formula='H', name='H+', compartment='Met_2c', charge=1)
pi_Met_2c = Metabolite('pi_Met_2c', formula='HO4P', name='xylose-D', compartment='Met_2c', charge=-2)
h2o_Met_2c = Metabolite('h2o_Met_2c', formula='H2O', name='H2O', compartment='Met_2c', charge=0)
actp_Met_2c = Metabolite('actp_Met_2c', formula='C2H3O5P', name='Acetyl phosphate', compartment='Met_2c', charge=-2)
for_Met_2c = Metabolite('for_Met_2c', formula='CHO2', name='Formate', compartment='Met_2c', charge= -1)
accoa_Met_2c = Metabolite('accoa_Met_2c', formula='C23H34N7O17P3S', name='Acetyl-CoA', compartment='Met_2c', charge=-4)
coa_Met_2c = Metabolite('coa_Met_2c', formula='C21H32N7O16P3S', name='Coenzyme A', compartment='Met_2c', charge=-4)
co_Met_2c = Metabolite('co_Met_2c', formula='CO', name='Carbon monoxide', compartment='Met_2c', charge=0)
fdred_Met_2c = Metabolite('fdred_Met_2c', formula='Fe8S8X', name='Ferredoxin (reduced) 2[4Fe-4S]', compartment='Met_2c', charge= -2)
fdox_Met_2c = Metabolite('fdox_Met_2c', formula='Fe8S8X', name='Ferredoxin (oxidized) 2[4Fe-4S]', compartment='Met_2c', charge= 0)
co2_Met_2c = Metabolite('co2_Met_2c', formula='CO2', name='CO2', compartment='Met_2c', charge= 0)
cfesp_Met_2c = Metabolite('cfesp_Met_2c', formula='C19CoN4R21', name='Corrinoid Iron sulfur protein', compartment='Met_2c', charge=-1)
mecfsp_Met_2c = Metabolite('mecfsp_Met_2c', formula='C20H3CoN4R21', name='Methylcorrinoid iron sulfur protein', compartment='Met_2c', charge=0)
_5mthf_Met_2c = Metabolite('_5mthf_Met_2c', formula='C20H24N7O6', name='5-Methyltetrahydrofolate', compartment='Met_2c', charge=-1)
mlthf_Met_2c = Metabolite('mlthf_Met_2c', formula='C20H21N7O6', name='5,10-Methylenetetrahydrofolate', compartment='Met_2c', charge=-2)
methf_Met_2c = Metabolite('methf_Met_2c', formula='C20H20N7O6', name='5,10-Methenyltetrahydrofolate', compartment='Met_2c', charge=-1)
thf_Met_2c = Metabolite('thf_Met_2c', formula='C19H21N7O6', name='5,6,7,8-Tetrahydrofolate', compartment='Met_2c', charge=-2)
_10fthf_Met_2c = Metabolite('_10fthf_Met_2c', formula='C20H21N7O7', name='10-Formyltetrahydrofolate', compartment='Met_2c', charge=-2)
h2_Met_2c = Metabolite('h2_Met_2c', formula='H2', name='Hydrogen', compartment='Met_2c', charge= 0)
nad_Met_2c = Metabolite('nad_Met_2c', formula='C21H26N7O14P2', name='Nicotinamide adenine dinucleotide', compartment='Met_2c', charge=-1)
nadh_Met_2c = Metabolite('nadh_Met_2c', formula='C21H27N7O14P2', name='Nicotinamide adenine dinucleotide - reduced', compartment='Met_2c', charge=-2)
h_Met_2i = Metabolite('h_Met_2i', formula='H', name='H+', compartment='Met_2i', charge=1)
fdxo_42_Met_2c = Metabolite('fdxo_42_Met_2c', formula='Fe12S12X', name='Ferredoxin - oxidized', compartment='Met_2c',charge=0)
fdxr_42_Met_2c = Metabolite('fdxr_42_Met_2c', formula='Fe12S12X', name='Ferredoxin - reduced', compartment='Met_2c',charge=-3)

mfr_b_Met_2c = Metabolite('mfr_b_Met_2c', formula='C34H44N6O14R', name='Methanofuran b', compartment='Met_2c', charge=-1)
formmfr_b_Met_2c = Metabolite('formmfr_b_Met_2c', formula='C35H43N6O15R', name='Formylmethanofuran b', compartment='Met_2c',charge=-2)
formh4spt_Met_2c = Metabolite('formh4spt_Met_2c', formula='C36H48N7O20P1', name='Formyltetrahydrosarcinapterin', compartment='Met_2c', charge=0)
h4spt_Met_2c = Metabolite('h4spt_Met_2c', formula='C35H48N7O19P1', name='Tetrahydrosarcinapterin', compartment='Met_2c',charge=0)
menylh4spt_Met_2c = Metabolite('menylh4spt_Met_2c', formula='C36H47N7O19P1', name='Methenyl-tetrahydrosarcinapterin',compartment='Met_2c', charge=1)
f420_2_Met_2c = Metabolite('f420_2_Met_2c', name='Coenzyme ferredoxin 420-2 (oxidized)', formula='C29H34N5O18P',compartment='Met_2c', charge=0)
f420_2h2_Met_2c = Metabolite('f420_2h2_Met_2c', name='Coenzyme ferredoxin 420-2 (reduced)', formula='C29H36O18N5P1',compartment='Met_2c', charge=0)
mleneh4spt_Met_2c = Metabolite('mleneh4spt_Met_2c', name='N5,N10-methylee-5,6,7,8-tetrahydromethanopterin', formula='C36H48N7O19P1', compartment='Met_2c', charge=0)
com_Met_2c = Metabolite('com_Met_2c', name='Coenzyme m', formula='C2H5O3S2', compartment='Met_2c', charge=0)
mcom_Met_2c = Metabolite('mcom_Met_2c', name='Methylcoenzyme m', formula='C3O3S2H7', compartment='Met_2c', charge=0)
mphen_Met_2c = Metabolite('mphen_Met_2c', name='Methanophenazine (oxidized)', formula='C37N2O1H50', compartment='Met_2c',charge=0)
mphenh2_Met_2c = Metabolite('mphenh2_Met_2c', name='Methanophenazine (reduced)', formula='C37N2O1H52', compartment='Met_2c',charge=0)
mh4spt_Met_2c = Metabolite('mh4spt_Met_2', name='N5-methyl-tetrahydrosarcinapterin', formula='C36H50N7O19P1',compartment='Met_2c', charge=0)

ch4_Met_2c = Metabolite('ch4_Met_2c', formula='CH4', name='Methane', compartment='Met_2c', charge=0)
cob_Met_2c = Metabolite('cob_Met_2c', formula='C11H19N1O7P1S1', name='Conenzyme B', compartment='Met_2c', charge=0)
hsfd_Met_2c = Metabolite('hsfd_Met_2c', formula='C13H22N1O10P1S3', name='Heterodisulfide', compartment='Met_2c', charge=0)

#FdhAB
reaction = Reaction('Met_2_FdhAB')  # What do we do with R's if they're not balancing?
reaction.name = 'Met_2: FdhAB'
reaction.subsystem = 'Met_2'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({for_Met_2c: -1.0,
                          h_Met_2c: -1.0,
                          f420_2_Met_2c: -1.0,
                          co2_Met_2c: 1.0,
                          f420_2h2_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# [1] Formylmethanofuran Dehydrogenase - FMFD_b
reaction = Reaction('Met_2_FwdA-H')  # What do we do with R's if they're not balancing?
reaction.name = 'Met_2: Formylmethanofuran Dehydrogenase'
reaction.subsystem = 'Met_2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({co2_Met_2c: -1.0,
                          fdred_Met_2c: -1.0,
                          h_Met_2c: -1.0,
                          mfr_b_Met_2c: -1.0,
                          fdox_Met_2c: 1.0,
                          h2o_Met_2c: 1.0,
                          formmfr_b_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# [2] Formylmethanofuran-tetrahydromethanopterin N-formyltransferase - FMFTSPFT_b
reaction = Reaction('Met_2_Ftr')
reaction.name = 'Met_2: Formylmethanofuran-Tetrahydromethanopterin N-Formyltransferase'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_2c: -1.0,
                          h4spt_Met_2c: -1.0,
                          mfr_b_Met_2c: 1.0,
                          formmfr_b_Met_2c: -1.0,
                          formh4spt_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# [3]? MTSPC
reaction = Reaction('Met_2_Mch')
reaction.name = 'Met_2: Mch'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_2c: -1.0,
                          h2o_Met_2c: 1.0,
                          menylh4spt_Met_2c: 1.0,
                          formh4spt_Met_2c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4RHi [9]
reaction = Reaction('Met_2_FrhABG')
reaction.name = 'Met_2: FrhABG'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({h2_Met_2c: -1.0,
                          f420_2_Met_2c: -1.0,
                          f420_2h2_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4MTSPD [4]
reaction = Reaction('Met_2_Mtd')
reaction.name = 'Met_2: Mtd'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_2c: 1.0,
                          f420_2_Met_2c: 1.0,
                          f420_2h2_Met_2c: -1.0,
                          mleneh4spt_Met_2c: 1.0,
                          menylh4spt_Met_2c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# F4MTSPR [5]
reaction = Reaction('Met_2_Mer')
reaction.name = 'Met_2: Mer'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({f420_2_Met_2c: 1.0,
                          mh4spt_Met_2c: 1.0,
                          f420_2h2_Met_2c: -1.0,
                          mleneh4spt_Met_2c: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# MTSPCMMT [6]
reaction = Reaction('Met_2_MtrA-H')
reaction.name = 'Met_2: MtrA-H'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

reaction.add_metabolites({h_Met_2c: -2.0,
                          h_Met_2i: 2.0,
                          h4spt_Met_2c: 1.0,
                          mh4spt_Met_2c: -1.0,
                          com_Met_2c: -1.0,
                          mcom_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# MCR
reaction = Reaction('Met_2_McrABCDG')
reaction.name = 'Met_2: McrABCDG'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ch4_Met_2c: 1.0,
                          cob_Met_2c: -1.0,
                          mcom_Met_2c: -1.0,
                          hsfd_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# HDR
reaction = Reaction('Met_2_MvhADG+HdrABC')
reaction.name = 'Met_2: HDR'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({h2_Met_2c: -2.0,
                          h_Met_2c: 2.0,
                          hsfd_Met_2c: -1.0,
                          cob_Met_2c: 1.0,
                          com_Met_2c: 1.0,
                          fdox_Met_2c: -1.0,
                          fdred_Met_2c: 1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

###Hydrogen Production###

reaction = Reaction('Met_2_ECH')
#The reaction in BiGG uses a different ferredoxin
#BiGG reaction is not balanced for H
reaction.name = 'Met_2: Energy conserving hydrogenase'
reaction.subsystem = 'Hydrogen Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({fdred_Met_2c: -1.0,
                          h_Met_2c: -4.0,
                          h2_Met_2c: 1.0,
                          fdox_Met_2c: 1.0,
                          h_Met_2i: 2.0})

model.add_reactions([reaction])

print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_2_ATPS4r')
#This reaction differs from the BiGG reaction because this model assumes a different compartment for ion motive force generation
reaction.name = 'Met_2: ATP Synthase'
reaction.subsystem = 'Energy Generation'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({adp_Met_2c: -1.0,
                          pi_Met_2c: -1.0,
                          h_Met_2i: -4.0,
                          atp_Met_2c: 1.0,
                          h_Met_2c: 3.0,
                          h2o_Met_2c: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#ATP Hydrolysis

#atp_Met_2c + h2o_Met_2c <-> adp_Met_2c + pi_Met_2c + h_Met_2c + ATP_COMM_e

reaction = Reaction('Met_2_ATP_Hydrolysis')
reaction.name = 'Met_2: ATP Hydrolysis'
reaction.subsystem = 'ATP Hydrolysis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

reaction.add_metabolites({atp_Met_2c: -1.0,
                          h2o_Met_2c: -1.0,
                          adp_Met_2c: 1.0,
                          pi_Met_2c: 1.0,
                          h_Met_2c: 1.0,
                          atp_COMM: 1.0})

model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# Met_2 - TRANSPORT RXNS

reaction = Reaction('Met_2_EX_H2')
reaction.name = 'Met_2: EX - H2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2_Met_2c: 1.0,
                          h2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('Met_2_EX_Protons')
reaction.name = 'Met_2: Transport - Protons'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h_Met_2c: 1.0,
                          h_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('Met_2_EX_ch4')
reaction.name = 'Met_2: EX - ch4'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ch4_Met_2c: 1.0,
                          ch4_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

reaction = Reaction('Met_2_EX_h2o')
reaction.name = 'Met_2: EX - h2o'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2o_Met_2c: 1.0,
                          h2o_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('Met_2_EX_co2')
reaction.name = 'Met_2: EX - co2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({co2_Met_2c: 1.0,
                          co2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


reaction = Reaction('Met_2_EX_for')
reaction.name = 'Met_2: EX - for'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({for_Met_2c: 1.0,
                          h_Met_2c: 1.0,
                          h_e: -1.0,
                          for_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# ACETATE EXCHANGE
reaction = Reaction('EX_ac')
reaction.name = 'Exchange - Acetate'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ac_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# H EXCHANGE
reaction = Reaction('EX_h')
reaction.name = 'Exchange - h'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

#H2 EXCHANGE
reaction = Reaction('EX_H2')
reaction.name = 'Exchange - H2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# FOR EXCHANGE
reaction = Reaction('EX_for')
reaction.name = 'Exchange - for'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({for_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# CO2 EXCHANGE
reaction = Reaction('EX_co2')
reaction.name = 'Exchange - co2'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({co2_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# H2O EXCHANGE
reaction = Reaction('EX_h2o')
reaction.name = 'Exchange - h2o'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({h2o_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# CH4 EXCHANGE
reaction = Reaction('EX_ch4')
reaction.name = 'Exchange - ch4'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ch4_e: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))

# ATP COMM EXCHANGE
reaction = Reaction('EX_ATP_COMM')
reaction.name = 'Exchange - Community ATP'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({atp_COMM: -1.0})
model.add_reactions([reaction])
print(reaction.name + ": " + str(reaction.check_mass_balance()))


Constraint_ATP_1 = model.problem.Constraint(model.reactions.DTU068_1_ATP_Hydrolysis.flux_expression  - 1.86513*model.reactions.Met_2_ATP_Hydrolysis.flux_expression, lb=0, ub=0)
model.add_cons_vars(Constraint_ATP_1)

Constraint_ATP_1 = model.problem.Constraint(model.reactions.Met_1_ATP_Hydrolysis.flux_expression  - 1.889665*model.reactions.DTU068_1_ATP_Hydrolysis.flux_expression, lb=0, ub=0)
model.add_cons_vars(Constraint_ATP_1)

#######################################################################################################################
########################################################################################################################

medium = model.medium

medium["EX_H2"] = 0
medium["EX_for"] = 0
medium["EX_co2"] = 0
medium["EX_ac"] = 1

model.medium = medium
print(model.medium)

# Add additional constraints

model.reactions.DTU068_1_ECH.lower_bound = 0
model.reactions.DTU068_1_HYDABC.lower_bound = 0
model.reactions.DTU068_1_NiFeH2ASE.lower_bound = 0

model.objective = 'EX_ATP_COMM'

pfba_solution = cobra.flux_analysis.pfba(model)
model.summary()
print(pfba_solution.fluxes)

pfba_solution.fluxes.to_json("fluxes.json")

workbook = xlsxwriter.Workbook('FluxResults.xlsx')
worksheet = workbook.add_worksheet('FluxResults')

worksheet.write(0, 0, "reaction")
worksheet.write(0, 1, "flux")

row = 1
col = 0

for rxn in model.reactions:
    worksheet.write(row, col, str(rxn))
    row += 1

row = 1
col += 1

for flux in pfba_solution.fluxes:
    worksheet.write(row, col, flux)
    row += 1

workbook.close()

cobra.io.save_json_model(model, "model.json")

print(f'{len(model.reactions)} reactions')
print(f'{len(model.metabolites)} metabolites')

cobra.io.write_sbml_model(model, "model.xml")


