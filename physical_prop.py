#from thermo import Mixture
#import fluids
#import math
import pandas as pd
import numpy as np
import streamlit as st
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashPureVLS,IAPWS95Gas,IAPWS95Liquid
from thermo.interaction_parameters import IPDB
#from thermo.property_package import GceosBase

gases_list = ['water', 'hydrogen', 'nitrogen', 'carbon dioxide', 'hydrogen sulfide','methane',
'ethane', 'propane', 'isobutane', 'n-butane', 'isopentane', 'n-pentane', 'hexane',
'heptane', 'octane', 'nonane']
liquid_list= ['water', 'ethanol', 'methanol', 'acetic acid', 'propylene glycol', 'glycerol', 'dimethyl sulfoxide', 'benzene', 'toluene', 'xylene', 'acetone', 'butanol', 'pentanol', 'hexanol', 'heptanol', 'octanol', 'nonanol', 'decanol', 'ethylene glycol', 'diethylene glycol', 'propylene carbonate', 'tetrahydrofuran', 'acetonitrile', 'formamide', 'isopropyl alcohol', 'methyl ethyl ketone', 'dioxane', 'pyridine', 'hexamethylphosphoramide','dimethylamine', 'diethylamine', 'triethylamine', 'trimethylamine', 'ethanolamine', 'diethanolamine', 'triethanolamine', 'methyldiethanolamine', 'piperazine']
c = gases_list + liquid_list
c.remove('water')
if 'constants' not in st.session_state:
    st.session_state.constants, st.session_state.correlations = ChemicalConstantsPackage.from_IDs(c)
    
constants = st.session_state.constants
correlations = st.session_state.correlations
kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')

zs = [1/(len(c)) for i in range(len(c))]
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=correlations.HeatCapacityGases)
flasher = FlashPureVLS(constants, correlations, liquids=[liquid], gas=gas, solids=[])
T1 = 273.15+30
state_1 = flasher.flash(P=100000, T=T1,zs=zs)
    
properties = {    "Mass": "kg",    "Length": "m",    "Time": "s",    "Temperature": "°C",    "Heat capacity": "Kcal/(kg*°C)",    "Enthalpy": "KCal/kg",    "thermal conductivity": "W/(m*°C)",    "Mass flow rate": "kg/hr",    "viscosity": "cP",    "density": "kg/m³","Cv": "Kcal/(kg*°C)","Cp": "Kcal/(kg*°C)"}
s = pd.Series(properties)    

def density(sg,temperature):
                api = (141.5/sg) - 131.5
                if api <=14.9:
                    eta= 0.00035
                elif api <= 34.9:
                    eta= 0.0004
                elif api <= 50.9:
                    eta= 0.0005
                elif api <= 63.9:
                    eta= 0.0006
                elif api <= 78.9:
                    eta= 0.0007
                elif api <= 88.9:
                    eta= 0.0008
                elif api <= 93.9:
                    eta= 0.00035 
                elif api <= 100:
                    eta= 0.0009   
                else: st.write("please enter a valid Specific gravity or use user-defined compositions in Gases' option") 
                density = sg*1000*(1-eta*((temperature*1.8+32)-60))
                return density

def thermo_prop(sg,t,prop_calc_table):
        t = 1.8*t+32
        thermal_coductivity = (0.813/sg)*(1-0.0003*(t-32)) *0.14422790000000002
        cp = (1/np.sqrt(sg))*(0.388+0.00045*t)
        cv = cp-(0.09/sg)
        latent_heat = (1/sg)*(110.9-0.09*t) * 0.555927
        if prop_calc_table.loc['thermal conductivity','Method'] != 'Two Linear points' :
            prop_calc_table.loc['thermal conductivity','Calculated_properties']= thermal_coductivity
            prop_calc_table.loc['thermal conductivity','Method']= 'Bureau Report 1929'
        if prop_calc_table.loc['Cp','Method']!= 'Two Linear points' :
            prop_calc_table.loc['Cp','Calculated_properties']= cp
            prop_calc_table.loc['Cp','Method']= 'Bureau Report 1929'
        prop_calc_table.loc['Cv','Calculated_properties'] = cv
        prop_calc_table.loc['latent heat','Calculated_properties']= latent_heat
        prop_calc_table.loc[['Cv','latent heat'],'Method']= 'Bureau Report 1929'
        return prop_calc_table
def vis_1point(t,analysis_temp,analysis_mu,sg,unit):
    if not unit:
        T =t
        mu = analysis_mu
        c = -0.8696
        b= np.log10(mu)-c
        s = 0.2008*b+1.6180
        log_mu = (b/(1+((T-analysis_temp)/310.93))**s)+c
        mu_calc =10**log_mu
        mu_calc = mu_calc/(0.001*density(sg,T))
    else:
        T =t
        mu = analysis_mu/(0.001*density(sg,analysis_temp))
        c = -0.8696
        b= np.log10(mu)-c
        s = 0.2008*b+1.6180
        log_mu = (b/(1+((T-analysis_temp)/310.93))**s)+c
        mu_calc =(10**log_mu)*(0.001*density(sg,analysis_temp))
    return mu_calc
def thermo_prop_LorGas(type):
        props = ['Phase','Vapor Fraction','density','Molecular Weight', 'Cp','Cv','K (Cp/Cv)', 'thermal conductivity','viscosity','Compressibility factor']
        prop_calc_table = pd.DataFrame(index=props,columns=['Calculated_properties'])
        if type == 'Gas':
            try:
                # Define the pipe and conditions
                pressure = float(st.number_input('Pressure in kg/cm2.a'))*98066.5
                temperature_K = float(st.number_input('Temperature in C')) + 273.15
                composition = st.multiselect('Components', gases_list)
                composition_table = pd.DataFrame(index=composition,columns=['mole fraction%'])
                
                comp_table = st.experimental_data_editor(composition_table)
                mole_fractions = {comp_table.index[i]: comp_table['mole fraction%'].astype('float64')[i]/100 for i in range(len(comp_table.index))}
                if sum(comp_table['mole fraction%'].astype('float64')) == 100:
                        st.success('Composition in Mol. percent completed!', icon="✅")
                if st.button("Calculate", key = 'calculations_tablegas'):
                    if sum(comp_table['mole fraction%'].astype('float64')) == 100:
                        
                        zs = [mole_fractions[i] if i in mole_fractions.keys() else 0 for i in c]
                        
                        gas_mixture = flasher.flash(P=pressure, T=temperature_K,zs=zs)
                        if 'water' in mole_fractions.keys() and mole_fractions['water'] == 1 :
                            gas = IAPWS95Gas(T=temperature_K, P=pressure, zs=zs)
                            liq = IAPWS95Liquid(T=temperature_K, P=pressure, zs=zs)
                            flasher_new= FlashPureVLS(constants, properties, liquids=[liq], gas=gas, solids=[])
                            mix2 = flasher_new.flash(T=temperature_K, P=pressure, zs=zs)
                            gas_mixture = mix2 
                        
                        prop_calc_table.loc['Phase','Calculated_properties'] = gas_mixture.phase
                        prop_calc_table.loc['Vapor Fraction','Calculated_properties'] = gas_mixture.VF
                        prop_calc_table.loc['thermal conductivity','Calculated_properties'] = gas_mixture.k()
                        prop_calc_table.loc['density','Calculated_properties'] = gas_mixture.rho_mass()
                        prop_calc_table.loc['Cp','Calculated_properties'] = gas_mixture.Cp_mass()/4184
                        prop_calc_table.loc['Cv','Calculated_properties'] = gas_mixture.Cv_mass()/4184
                        prop_calc_table.loc['viscosity','Calculated_properties'] = gas_mixture.mu()*1000
                        prop_calc_table.loc['Molecular Weight','Calculated_properties'] = gas_mixture.MW()
                        prop_calc_table.loc['Compressibility factor','Calculated_properties'] = gas_mixture.Z()
                        prop_calc_table.loc['K (Cp/Cv)','Calculated_properties'] = gas_mixture.isentropic_exponent()
                        prop_calc_table.loc['Enthalpy','Calculated_properties'] = gas_mixture.H_mass()/4184
                        prop_calc_table.loc['LHV','Calculated_properties'] = -sum(pd.Series(zs)*pd.Series(gas_mixture.Hcs_lower_mass).fillna(0)*(pd.Series(gas_mixture.MWs)/gas_mixture.MW()))/4184
                        #st.write(-gas_mixture.Hc_lower_mass()) #/4184)
                        prop_calc_table = prop_calc_table.merge(s.rename('Units'), left_index=True,right_index=True, how='left')
                        prop_calc_table.loc[:,'Method']= 'Thermo Library'
                        
                        st.write(prop_calc_table)
                        
                            
                        
            except IndexError: pass
            except (ValueError): st.write('Please check your input')
             
        if type == 'Liquid':
            try:
                # Define the pipe and conditions
                pressure = float(st.number_input('Pressure in kg/cm2.a'))*98066.5
                temperature_K = float(st.number_input('Temperature in C')) + 273.15
                composition = st.multiselect('Components', liquid_list)
                composition_table = pd.DataFrame(index=composition,columns=['Volume fraction%'])
                
                comp_table = st.experimental_data_editor(composition_table)
                
                
                mole_fractions = {comp_table.index[i]: comp_table['Volume fraction%'].astype('float64')[i]/100 for i in range(len(comp_table.index))}
                if sum(comp_table['Volume fraction%'].astype('float64')) == 100:
                        st.success('Composition in Mol. percent completed!', icon="✅")
                if st.button("Calculate", key = 'calculations_tableliquid'):
                    if sum(comp_table['Volume fraction%'].astype('float64')) == 100:
                        
                        zs = [mole_fractions[i] if i in mole_fractions.keys() else 0 for i in c]
                        
                        mixture = flasher.flash(P=pressure, T=temperature_K,zs=zs)
                        #mixture_trial = Mixture([i for i in mole_fractions.keys()], ws=[i for i in mole_fractions.values()], T=temperature_K, P=pressure, pkg= GceosBase)
                        #st.write(mixture_trial.rho)
                        if 'water' in mole_fractions.keys() and mole_fractions['water'] == 1 :
                            gas = IAPWS95Gas(T=temperature_K, P=pressure, zs=zs)
                            liq = IAPWS95Liquid(T=temperature_K, P=pressure, zs=zs)
                            flasher_new= FlashPureVLS(constants, properties, liquids=[liq], gas=gas, solids=[])
                            mix2 = flasher_new.flash(T=temperature_K, P=pressure, zs=zs)
                            mixture = mix2 
                        if 'water' in mole_fractions.keys() and mole_fractions['water'] < 1 :
                            prop_calc_table.loc['density','Calculated_properties'] = mixture.rho_mass()*(1+mole_fractions['water']*(998/847.38))
                        prop_calc_table.loc['Phase','Calculated_properties'] = mixture.phase
                        prop_calc_table.loc['Vapor Fraction','Calculated_properties'] = mixture.VF
                        prop_calc_table.loc['thermal conductivity','Calculated_properties'] = mixture.k()
                        
                        if 'water' in mole_fractions.keys() and mole_fractions['water'] < 1 :
                            prop_calc_table.loc['density','Calculated_properties'] = mixture.rho_mass()*(1+mole_fractions['water']*((998-847.38)/847.38))
                            
                        prop_calc_table.loc['Cp','Calculated_properties'] = mixture.Cp_mass()/4184
                        prop_calc_table.loc['Cv','Calculated_properties'] = mixture.Cv_mass()/4184
                        prop_calc_table.loc['viscosity','Calculated_properties'] = mixture.mu()*1000
                        prop_calc_table.loc['Molecular Weight','Calculated_properties'] = mixture.MW()
                        prop_calc_table.loc['Compressibility factor','Calculated_properties'] = mixture.Z()
                        prop_calc_table.loc['K (Cp/Cv)','Calculated_properties'] = mixture.isentropic_exponent()
                        prop_calc_table.loc['Enthalpy','Calculated_properties'] = mixture.H_mass()/4184
                        prop_calc_table.loc['LHV','Calculated_properties'] = -sum(pd.Series(zs)*pd.Series(mixture.Hcs_lower_mass).fillna(0)*(pd.Series(mixture.MWs)/mixture.MW()))/4184
                        prop_calc_table = prop_calc_table.merge(s.rename('Units'), left_index=True,right_index=True, how='left')
                        prop_calc_table.loc[:,'Method']= 'Thermo Library'
                        
                        
                        #if 'water' in mole_fractions.keys():
                         #   mixture2 = Mixture([i for i in mole_fractions.keys()], ws=[i for i in mole_fractions.values()], T=temperature_K, P=pressure, pkg= GceosBase)
                          #  prop_calc_table.loc['density','Calculated_properties'] = mixture2.rho
                        st.write(prop_calc_table)

                        
            except IndexError: pass
            except (ValueError): st.write('Please check your input')
        
def main():
    
    phases  = st.selectbox('fluids presesnt',('Oil Fractions', 'liquid', 'Gas'), key='phases')
    if phases  == 'Gas':
        thermo_prop_LorGas('Gas')
    if phases == 'liquid':
        thermo_prop_LorGas('Liquid')
    elif phases  == 'Oil Fractions':
        
        try:
            props = ['Cp','Cv', 'thermal conductivity','latent heat','viscosity']
            prop_calc_table = pd.DataFrame(index=props,columns=['Calculated_properties','Method'])
            
            two_points = st.selectbox('Use 2 points of a certain property?',('No', 'Yes'), key='two_points')
            if two_points == "Yes":
                temperature = float(st.number_input('fluid Temperature in C', key='target_temp')) 
                prop_menu = st.multiselect('select property with two data points',['viscosity', 'specific gravity', 'Cp', 'thermal conductivity'])
                temperature1 = float(st.number_input('point 1 Temperature in C', key='T1')) 
                temperature2 = float(st.number_input('point 1 Temperature in C', key='T2')) 
                prop_table = pd.DataFrame(index=prop_menu,columns=['point 1','point 2'])
                prop_table_st = st.experimental_data_editor(prop_table)
                sg = float(st.number_input('Specific gravity at 15.56 C'))
                
            if two_points != 'Yes':           
                temperature = float(st.number_input('fluid Temperature in C', key='target_temp1'))
                sg = float(st.number_input('Specific gravity at 15.56 C'))
                vis_1point_select  = st.selectbox('Calculate viscosity using 1 point?',('No', 'Yes'), key='vis_1pointer')
                if vis_1point_select == 'Yes':
                    temperature_analysis = float(st.number_input('analysis Temperature in C', key='analysis_temp1'))
                    vis_analysis = float(st.number_input('Viscosity at analysis temperature in C.st', key='analysis_vis'))
                    unit = st.checkbox('viscosity Unit is in cP not C.st')
                    viscosity_calc = vis_1point(temperature,temperature_analysis,vis_analysis,sg,unit)
                    prop_calc_table.loc['viscosity','Calculated_properties'] = viscosity_calc
                    prop_calc_table.loc['viscosity','Method']= 'One point - A. Miadonye and V.R. Puttagunta'
                    
            
            prop_calc_table = thermo_prop(sg,temperature,prop_calc_table)
            prop_calc_table.loc['density'] = [density(sg,temperature),'Nelson']
        except (ZeroDivisionError,UnboundLocalError): pass
        if st.button("Calculate", key = 'calculations_tableLiq'):
                try:
                    if two_points == 'Yes':
                        for i in prop_menu:
                            
                            if i != 'viscosity':
                                
                                A = np.array([[1,temperature1], [1,temperature2]])
                                B = np.array([float(prop_table_st.loc[i,'point 1']),float(prop_table_st.loc[i,'point 2'])])
                                C = np.linalg.solve(A, B)
                                prop_calc_table.loc[i,'Calculated_properties']=C[0]+temperature*C[1]
                                prop_calc_table.loc[i,'Method']= 'Two Linear points'
                                
                                
                            else:
                    
                                # define the points (x1, y1) and (x2, y2)
                                x1 = temperature1+273.15
                                y1 = float(prop_table_st.loc[i,'point 1'])
                                x2 = temperature2+273.15
                                y2 = float(prop_table_st.loc[i,'point 2'])

                                # compute the values of z1 and z2
                                z1 = np.log10(y1)
                                z2 = np.log10(y2)

                                # solve for a and b
                                a = (z2 - z1) / (x2 - x1)
                                b = z1 - a * x1
                                viscosity = 10**(a*(temperature+273.15)+b)
                                # print the values of a and b
                                prop_calc_table.loc[i,'Calculated_properties'] = viscosity
                                prop_calc_table.loc[i,'Method']= 'Two Log points'
                                prop_calc_table = prop_calc_table.merge(s.rename('Units'), left_index=True,right_index=True).reindex(columns=['Calculated_properties', 'Units', 'Method'])
                        st.write(prop_calc_table.dropna(how='any'))
                    else:
                        
                        prop_calc_table = prop_calc_table.merge(s.rename('Units'), left_index=True,right_index=True).reindex(columns=['Calculated_properties', 'Units', 'Method'])
                        st.write(prop_calc_table.dropna(how='any'))
                except (ValueError,np.linalg.LinAlgError): st.write('Please check your points input')
if __name__ == '__main__':
    main()
 