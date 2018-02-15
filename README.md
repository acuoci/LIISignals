# LIISignals
Simulation of Laser-Induced Incandescence (LII) signals

LIISignal is a numerical library written in C++ for modeling the temporal behavior of LII signals by solving a coupled energy 
and mass balance including the absorption of laser energy, the heat loss due to vaporization from the particle surface, 
heat conduction to the surrounding gas, and heat loss due to radiation, according to the model developed by Hofman et al. (2007).

A detailed description of governing equations is available in the following paper:

**M. Hofmann,  B. Kock, C. Schulz**  
*A web-based interface for modeling laser-induced incandescence (LIISim)*  
Proceedings of the European Combustion Meeting (2007)  
http://web.liisim.com/Hofmann_et_al_2007_LIISim.pdf

The present implementation of the model tries to reproduce the web-based interface developed by Hofmann et al. 
and available at the following web address:
http://web.liisim.com/
