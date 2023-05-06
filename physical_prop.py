from thermo import Mixture
import fluids
import math
import pandas as pd
import numpy as np
import streamlit as st

gases_list = ['water', 'hydrogen', 'nitrogen', 'carbon dioxide', 'hydrogen sulfide','methane',
'ethane', 'propane', 'isobutane', 'n-butane', 'isopentane', 'n-pentane', 'hexane',
'heptane', 'octane', 'nonane']
liquid_list= ['water', 'ethanol', 'methanol', 'acetic acid', 'propylene glycol', 'glycerol', 'dimethyl sulfoxide', 'benzene', 'toluene', 'xylene', 'acetone', 'butanol', 'pentanol', 'hexanol', 'heptanol', 'octanol', 'nonanol', 'decanol', 'ethylene glycol', 'diethylene glycol', 'propylene carbonate', 'tetrahydrofuran', 'acetonitrile', 'formamide', 'isopropyl alcohol', 'methyl ethyl ketone', 'dioxane', 'pyridine', 'hexamethylphosphoramide','dimethylamine', 'diethylamine', 'triethylamine', 'trimethylamine', 'ethanolamine', 'diethanolamine', 'triethanolamine', 'methyldiethanolamine', 'piperazine']
c = gases_list + liquid_list
mixture = Mixture(c, zs=[1/(len(c)) for i in range(len(c))], T=298, P=100000)
#['water','ethanol','dimethylamine', 'diethylamine', 'triethylamine', 'trimethylamine', 'ethanolamine', 'diethanolamine', 'triethanolamine', 'methyldiethanolamine', 'piperazine']
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
def vis_1point(t,analysis_temp,analysis_mu):
    T =t
    mu = analysis_mu
    c = -0.8696
    b= np.log10(mu)-c
    s = 0.2008*b+1.6180
    log_mu = (b/(1+((T-analysis_temp)/310.93))**s)+c
    mu_calc =10**log_mu
    return mu_calc
def thermo_prop_LorGas(type):
        props = ['density', 'Cp','Cv', 'thermal conductivity','latent heat','viscosity']
        prop_calc_table = pd.DataFrame(index=props,columns=['Calculated_properties','Method'])
        if type == 'Gas':
            try:
                # Define the pipe and conditions
                pressure = float(st.number_input('Pressure in kg/cm2.a'))*98066.5
                temperature_K = float(st.number_input('Temperature in C')) + 273.15
                composition = st.multiselect('Components', gases_list)
                composition_table = pd.DataFrame(index=composition,columns=['mole fraction%'])
                
                comp_table = st.experimental_data_editor(composition_table)
                mole_fractions = {comp_table.index[i]: comp_table['mole fraction%'].astype('float64')[i]/100 for i in range(len(comp_table.index))}
                if st.button("Calculate", key = 'calculations_tablegas'):
                    if sum(comp_table['mole fraction%'].astype('float64')) == 100:
                        
                        #mole_fractions =  {"methane": 0.8, "ethane": 0.2}
                        gas_mixture = Mixture(list(mole_fractions.keys()), zs=list(mole_fractions.values()), T=temperature_K, P=pressure)
                        
                        if gas_mixture.phase == 'g':
                            prop_calc_table.loc['thermal conductivity','Calculated_properties'] = gas_mixture.kg
                            prop_calc_table.loc['density','Calculated_properties'] = gas_mixture.rho
                            prop_calc_table.loc['Cp','Calculated_properties'] = gas_mixture.Cp/4184
                            prop_calc_table.loc['Cv','Calculated_properties'] = gas_mixture.Cvg/4184
                            prop_calc_table.loc['viscosity','Calculated_properties'] = gas_mixture.mu*1000
                            prop_calc_table.loc['Molecular Weight','Calculated_properties'] = gas_mixture.MWg
                            prop_calc_table.loc['Compressibility factor','Calculated_properties'] = gas_mixture.Z
                            prop_calc_table.loc['K (Cp/Cv)','Calculated_properties'] = gas_mixture.isentropic_exponent
                            prop_calc_table.loc[:,'Method']= 'Thermo Library'
                            st.write(prop_calc_table)
                        else: st.warning('Liquid phase presence in fluid')
                        
            except IndexError: pass
        if type == 'Liquid':
            try:
                # Define the pipe and conditions
                pressure = float(st.number_input('Pressure in kg/cm2.a'))*98066.5
                temperature_K = float(st.number_input('Temperature in C')) + 273.15
                composition = st.multiselect('Components', liquid_list)
                composition_table = pd.DataFrame(index=composition,columns=['weight fraction%'])
                
                comp_table = st.experimental_data_editor(composition_table)
                
                
                mole_fractions = {comp_table.index[i]: comp_table['weight fraction%'].astype('float64')[i]/100 for i in range(len(comp_table.index))}
                if st.button("Calculate", key = 'calculations_tableliquid'):
                    if sum(comp_table['weight fraction%'].astype('float64')) == 100:
                        
                        
                        gas_mixture = Mixture(list(mole_fractions.keys()), ws=list(mole_fractions.values()), T=temperature_K, P=pressure)
                        
                        prop_calc_table.loc['thermal conductivity','Calculated_properties'] = gas_mixture.kl
                        prop_calc_table.loc['density','Calculated_properties'] = gas_mixture.rhol
                        prop_calc_table.loc['Cp','Calculated_properties'] = gas_mixture.Cpl/4184
                        prop_calc_table.loc['Cv','Calculated_properties'] = gas_mixture.Prl/4184
                        prop_calc_table.loc['viscosity','Calculated_properties'] = gas_mixture.mul*1000
                        #prop_calc_table.loc['latent heat','Calculated_properties'] = gas_mixture.Hvaps/4184
                        prop_calc_table.loc[:,'Method']= 'Thermo Library'
                        st.write(prop_calc_table)
            except IndexError: pass
def main():
    
    phases  = st.selectbox('fluids presesnt',('Liquid H.Cs', 'liquid', 'Gas'), key='phases')
    if phases  == 'Gas':
        thermo_prop_LorGas('Gas')
    if phases == 'liquid':
        thermo_prop_LorGas('Liquid')
    elif phases  == 'Liquid H.Cs':
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
                    vis_analysis = float(st.number_input('Viscosity at analysis temperature', key='analysis_vis'))
                    viscosity_calc = vis_1point(temperature,temperature_analysis,vis_analysis)
                    prop_calc_table.loc['viscosity','Calculated_properties'] = viscosity_calc
                    prop_calc_table.loc['viscosity','Method']= 'One point - A. Miadonye and V.R. Puttagunta'
                    
            
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
            else: st.write("please enter a valid Specific gravity") 
            density = sg*1000*(1-eta*((temperature*1.8+32)-60))
            prop_calc_table = thermo_prop(sg,temperature,prop_calc_table)
            prop_calc_table.loc['Density'] = [density,'Nelson']
        except ZeroDivisionError: pass
        if st.button("Calculate", key = 'calculations_tableLiq'):
                
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
                    st.write(prop_calc_table)
                else: st.write(prop_calc_table)
if __name__ == '__main__':
    main()
