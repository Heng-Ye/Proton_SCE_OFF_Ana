//Basic parameters ------------------------------------------//
//PDG
//const int pdg=211; //pi+
const int pdg=2212; //proton
const double m_proton=0.938272046; //proton_mass, unit:GeV/c2

//TPC boundary
double minX =  -360.0;
double maxX = 360.0;
double minY =0.0;
double maxY = 600.0;
double minZ =  0.0; // G10 at z=1.2
double maxZ = 695.0;

//some constants
double NA=6.02214076e23;
double MAr=39.95; //gmol
double Density = 1.39; // g/cm^3

//Misc. Cut values --------------------------//
//XY-Cut
//double mean_x=-29.73; //prod4 mc
//double mean_y=422.4;
//double dev_x=1.5*4.046;
//double dev_y=1.5*3.679;

//double mean_x_data=-26.58; //prod4 data
//double mean_y_data=423.5; //prod4 data
//double dev_x_data=1.5*3.753; //prod4 data
//double dev_y_data=1.5*4.354; //prod4 data

double mean_x=-29.25; //prod4a mc
double mean_y=422.1; //prod4a mc
double dev_x=1.5*4.4; //prod4a mc
double dev_y=1.5*3.883; //prod4a mc

double mean_x_data=-26.59; //prod4a data
double mean_y_data=423.5; //prod4a data
double dev_x_data=1.5*3.744; //prod4a data
double dev_y_data=1.5*4.364; //prod4a data


//dedx cut
double dedx_min=30.;

//beam quality cut -----------------------------------------------------------------//
double min1_dx=0.; //new p3
double min2_dx=2.0; //new p3
double min1_dy=-2.0; //new p3
double min2_dy=2.0; //new p3
double min1_dz=-.5; //new p3
double min2_dz=1.5; //new p3

double mean_StartZ=5.01422e-01; //prod4a sce off
double sigma_StartZ=7.02285e-02; //prod4a sce off
double mean_StartY=4.22214e+02; //prod4a sce off
double sigma_StartY=4.22176e+00; //prod4a sce off
double mean_StartX=-3.09019e+01; //prod4a sce off
double sigma_StartX=4.64495e+00; //prod4a sce off

//double min1_z=5.10816e-02-3.*2.13366e-01; //p4 
//double min2_z=5.10816e-02+3.*2.13366e-01; //p4
//double min1_y=4.21863e+02-3.*4.11359e+00; //p4
//double min2_y=4.21863e+02+3.*4.11359e+00; //p4
//double min1_x=-3.05895e+01-3.*4.69242e+00; //p4  
//double min2_x=-3.05895e+01+3.*4.69242e+00; //p4

double min1_z=mean_StartZ-3.*sigma_StartZ; //p4a sce off
double min2_z=mean_StartZ+3.*sigma_StartZ; //p4a sce off
double min1_y=mean_StartY-3.*sigma_StartY; //p4a sce off
double min2_y=mean_StartY+3.*sigma_StartY; //p4a sce off
double min1_x=mean_StartX-3.*sigma_StartX; //p4a sce off
double min2_x=mean_StartX+3.*sigma_StartX; //p4a sce off

double dx_min=-3.; double dx_max=3.;
double dy_min=-3.; double dy_max=3.;
double dz_min=-3.; double dz_max=3.;
//double dxy_min=-1.; double dxy_max=3.;
double dxy_min=-1.; double dxy_max=4.;

double costh_min = 0.96; double costh_max = 2; //prod4a

double cosine_beam_primtrk_min=(9.92194e-01)-4.*(3.96921e-03); //p4 [spec]


//stopping proton cut
//double mean_norm_trklen_csda=9.32064e-01; //prod4 spec
//double sigma_norm_trklen_csda=1.32800e-02; //prod4 spec
//double mean_norm_trklen_csda=9.01289e-01; //prod4a (spec)
//double sigma_norm_trklen_csda=7.11431e-02; //prod4a (spec)
double mean_norm_trklen_csda=8.90792e-01; //prod4a (spec)
double sigma_norm_trklen_csda=7.02755e-02; //prod4a (spec)
double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;

//Mean         8.90792e-01   1.58432e-03   9.74031e-07   9.17169e-02
//Sigma        7.02755e-02   1.31717e-03   3.22586e-06   3.08065e-02

//new inel cut using chi^2 info
double pid_1=10.;
double pid_2=10.;

//Mean         3.99025e+00   5.67304e-02   5.11613e-06  -2.13497e-05
//Sigma        1.39141e+00   5.97052e-02   6.70322e-06  -5.50969e-03
