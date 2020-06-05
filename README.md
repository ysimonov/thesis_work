# thesis_work
Programs for my thesis

The work is based on statistical characterization of multipath fading channel. 
Particularly, my research was directed towards understanding of statistical properties of 
the source of electromagnetic radiation. 
The scripts allow efficient processing of S-parameter data, calculation of system channel transfer function
under the assumption that source and load impedances are perfectly matched with 50 Ohm characteristic impedance
of the transmission line.
The channel sounding was based on frequency sweep of 1-9GHz UWB channel formed by multipath medium, and two
BBHA 9120 D Double Ridge Broadband Horn antennas. FieldFox N9915A 9GHz Microwave Analyzer was used for S-parameter test.
Transfer function is fitted by a rational fit and is used to estimated impulse response, assuming the system is LTI. 
The order of rational polynomial was chosen > 200 Degree to achieve sufficient accuracy. 

I also wrote a separate subroutine that performs sampled continuous inverse Fourier transform and Fourier transform of 
recorded data, which closely replicates time-domain measurement of the Vector Network Analyzer.
The impulse response was convolved with thermal, AM and PM signals to understand how the radio waves will change when
subjected to multipath medium. 
