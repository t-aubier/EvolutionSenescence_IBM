last update: June 23th, 2022


#### General Information ####

1. Title of Code: Code from: Deleterious mutation accumulation and the evolution of senescence

2. Author Information
	Corresponding Investigator
		Name: Dr Thomas G. Aubier
		Institution: University of North Carolina at Chapel Hill, USA
		Email: thomas.aubier@normalesup.org

	Co-investigator 1
		Name: Dr Matthias Galipaud
		Institution: University of Zürich, Switzerland


3. Date of code development: 2020-2022

4. Funding sources that supported code development: Swiss National Science Foundation, and National Science Foundation

5. Recommended citation for this code: Aubier and Galipaud (2022), Code from: Deleterious mutation accumulation and the evolution of senescence, GitHub repository, https://github.com/t-aubier/EvolutionSenescence_IBM



#### System Requirements ####

Simulations are run using C++ (compiler g++ 7.5.0).



#### Installation guide ####

Install C++: https://code.visualstudio.com/docs/languages/cpp

Also install the Boost C++ Libraries development files.
On linux: sudo apt-get install libboost-all-dev

Overall installation time is short: < 5 min



#### Script description ####

Oriented-based programming in C++ with modeling of two different classes 'Population' and 'Individual':

"ClassPopulation.cpp" and "ClassPopulation.h": script defining the class 'Population' referring to a population of individuals; see the header files "ClassPopulation.h" for descriptions of the variables and the functions.

"ClassIndividual.cpp" and "ClassIndividual.h": script defining the class 'ClassIndividual' referring a single individual; see the header files "ClassIndividual.h" for descriptions of the variables and the functions.

"OtherFunctions.cpp" and "OtherFunctions.h": script defining additional basic functions (e.g., calculation of mean, generation of random numbers); see the header files "OtherFunctions.h" for descriptions of the functions.

The script "Simulation_Accumulation.cpp" is used to run forward simulations. All variables are described in the functions 'main' within this file.

The script "Sensitivity_Accumulation.cpp" is used to run a sensitivity analysis (variations of mortality and birth rates). All variables are described in the functions 'main' within this file.

In the folder "AnalysisR", the R script to draw graphs resulting from the simulations using R language.

TimeSeries.R                    Plot figures describing the time series; obtained using 'SimulationNetwork.cpp'
SensitivityAnalysis.R           Plot figures describing the outcome of multiple simulations (like in Fig. 4); obtained using 'SimulationNetwork.cpp'
AnalysisMutant.R                Plot figures describing the fate of mutants with specific trait values (like in Fig. 5); obtained using 'AnalysisMutant.cpp'



#### Demo and instructions for use ####

All functions are written in the script files "ClassIndividual.cpp" and "ClassPopulation.cpp" (see the header files "ClassIndividual.h" and "ClassPopulation.h" for descriptions of the variables and the functions). This file must remain untouched. Additional functions are written in the files "OtherFunctions.cpp" and "OtherFunctions.h". All these files must remain untouched.

Simulations are run using the scripts "Simulation_Accumulation.cpp" or "Sensitivity_Accumulation.cpp". One can change parameter values there. The command line (in Linux) is:

g++ -std=c++14 -O2  Simulation_Accumulation.cpp ClassPopulation.cpp ClassIndividual.cpp OtherFunctions.cpp
./a.out

g++ -std=c++14 -O2  Sensitivity_Accumulation.cpp ClassPopulation.cpp ClassIndividual.cpp OtherFunctions.cpp
./a.out

Data are then stored in the folder "Data".

In the folder "AnalysisR", use the R script to draw graphs resulting from the simulations using R language.

TimeSeries.R                    Plot figures describing the time series; obtained using 'SimulationNetwork.cpp'
SensitivityAnalysis.R           Plot figures describing the outcome of multiple simulations (like in Fig. 4); obtained using 'SimulationNetwork.cpp'
AnalysisMutant.R                Plot figures describing the fate of mutants with specific trait values (like in Fig. 5); obtained using 'AnalysisMutant.cpp'

Note that to draw the figures shown in the manuscript, considerable computing power was needed (via the use of a computing cluster; each simulation typically lasted ~24 hours). That is why the figures created by the Demo code are different from the ones in the manuscript. To reproduce the figures in the manuscript, one can change the parameter values to match the default values (default value are commented in the code).
