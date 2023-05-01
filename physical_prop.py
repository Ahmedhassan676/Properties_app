from thermo import Mixture
import fluids
import math
import pandas as pd
import numpy as np
import streamlit as st
import ht
gases_list = ['water', 'hydrogen', 'helium', 'nitrogen', 'carbon dioxide', 'hydrogen sulfide','methane',
'ethane', 'propane', 'isobutane', 'n-butane', 'isopentane', 'n-pentane', 'hexane',
'heptane', 'octane', 'nonane']
def thermo_prop(sg,t):
        t = 1.8*t+32
        thermo_dict = {}
        thermal_coductivity = (0.813/sg)*(1-0.0003*(t-32))
        cp = (1/np.sqrt(sg))*(0.388+0.00045*t)
        cv = cp-(0.09/sg)
        latent_heat = (1/sg)*(110.9-0.09*t)
        thermo_dict['thermal_conductivity'] = thermal_coductivity
        thermo_dict['Cp'] = cp
        thermo_dict['Cv'] = cv
        thermo_dict['latent_heat'] = latent_heat
        return thermo_dict
def main():
    phases = st.multiselect('fluids presesnt',['vapor H.Cs', 'Liquid H.Cs', 'Steam', 'water', 'Gases H.Cs'])
    if 'Gases H.Cs' in phases:
        try:
            # Define the pipe and conditions
            pressure = float(st.number_input('Pressure in kg/cm2.a'))*98066.5
            temperature_K = float(st.number_input('Temperature in C')) + 273.15
            composition = st.multiselect('gases composition', gases_list)
            composition_table = pd.DataFrame(index=composition,columns=['mole fraction%'])
            
            comp_table = st.experimental_data_editor(composition_table)
            
            st.write(comp_table['mole fraction%'][0])
            mole_fractions = {comp_table.index[i]: comp_table['mole fraction%'].astype('float64')[i]/100 for i in range(len(comp_table.index))}
            if sum(comp_table['mole fraction%'].astype('float64')) == 100:
                # Create a gas mixture
                st.write(mole_fractions)
                #mole_fractions =  {"methane": 0.8, "ethane": 0.2}
                gas_mixture = Mixture(list(mole_fractions.keys()), zs=list(mole_fractions.values()), T=temperature_K, P=pressure)
                thermo_dict_gas = {}
                thermo_dict_gas['thermal_conductivity'] = gas_mixture.kg
                thermo_dict_gas['density'] = gas_mixture.rho
                thermo_dict_gas['Cp'] = gas_mixture.Cp
                thermo_dict_gas['visocsity'] = gas_mixture.mu
                
                st.write(thermo_dict_gas)
        except IndexError: pass
    if 'Liquid H.Cs' in phases:
        temperature = float(st.number_input('Temperature in C')) 
        sg = float(st.number_input('Specific gravity'))
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
        st.write(thermo_prop(sg,temperature))
        st.write("density at {} equals {} kg/m3".format(temperature, density))
if __name__ == '__main__':
    main()
