#define ProtonOnlyTrackLength_cxx
#include "ProtonOnlyTrackLength.h"

#include <TH2.h>
#include <TH1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include "TVector3.h"

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

//#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
//#include "./headers/BasicFunctions.h"
#include "./headers/util.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
//#include "./headers/sce_map.h"

using namespace std;
using namespace ROOT::Math;


void ProtonOnlyTrackLength::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

	//book histograms --------------------------------------------------------------------------------------------//
	int n_b=150;
	//int n_b=30;
	double b_min=0;
	double b_max=150;

	TH1D *h1d_trklen_CaloSz=new TH1D(Form("h1d_trklen_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz=new TH1D(Form("h1d_endz_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_stz_CaloSz=new TH1D(Form("h1d_stz_CaloSz"), Form("MC CaloSz"), 500,-5,5);
	TH1D *h1d_sty_CaloSz=new TH1D(Form("h1d_sty_CaloSz"), Form("MC CaloSz"), 100,400,450);
	TH1D *h1d_stx_CaloSz=new TH1D(Form("h1d_stx_CaloSz"), Form("MC CaloSz"), 180,-60,30);
	TH1D *h1d_stxy_CaloSz=new TH1D(Form("h1d_stxy_CaloSz"), Form("MC CaloSz"), 500,-10,10);

	TH1D *h1d_trklen_Pos=new TH1D(Form("h1d_trklen_Pos"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_cosine_Pos=new TH1D(Form("h1d_cosine_Pos"), Form("MC CaloSz"), 100, .9, 1.);
	TH1D *h1d_trklen_BQ=new TH1D(Form("h1d_trklen_BQ"), Form("MC CaloSz"), n_b, b_min, b_max);
/*
	TH1D *h1d_trklen_CaloSz_el=new TH1D(Form("h1d_trklen_CaloSz_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_inel=new TH1D(Form("h1d_trklen_CaloSz_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midcosmic=new TH1D(Form("h1d_trklen_CaloSz_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midpi=new TH1D(Form("h1d_trklen_CaloSz_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midp=new TH1D(Form("h1d_trklen_CaloSz_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midmu=new TH1D(Form("h1d_trklen_CaloSz_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_mideg=new TH1D(Form("h1d_trklen_CaloSz_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midother=new TH1D(Form("h1d_trklen_CaloSz_midother"), Form("midother"), n_b, b_min, b_max);
*/

	//up-stream E-loss study
	float bx_min=-50;
	float bx_max=-10;
	int n_bx=40;
	float by_min=405;
	float by_max=440;
	int n_by=35;
	//TH2D *bx_by_RecoEl=new TH2D("bx_by_RecoEl","", n_bx, bx_min, bx_max, n_by, by_min, by_max);
	//TH2D *bx_by_RecoInEl=new TH2D("bx_by_RecoInEl","", n_bx, bx_min, bx_max, n_by, by_min, by_max);
	TH2D *bx_by_RecoAll=new TH2D("bx_by_RecoAll","", n_bx, bx_min, bx_max, n_by, by_min, by_max);
	//TH2D *bx_by_RecoMidP=new TH2D("bx_by_RecoMidP","", n_bx, bx_min, bx_max, n_by, by_min, by_max);

	TH1D *h1d_ntrklen_BQ=new TH1D(Form("h1d_ntrklen_BQ"), Form(""), 140,0,1.4);
	TH1D *h1d_pid_BQ=new TH1D(Form("h1d_pid_BQ"), Form(""), 500,0,100);


	//------------------------------------------------------------------------------------------------------------//
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//only select protons	
		if (beamtrackPdg!=pdg) continue; //only interested in protons

		//Event Selection Cut -- Part 1 ----------------------------------//
		bool IsBeamMatch=false; //if recostructed the right track (recoID=truthID)
		bool IsPandoraSlice=false; //pandora slice cut (can pandora reconstruct this track)
		bool IsCaloSize=false; //if calo size not empty
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (timeintersection->size()) IsIntersection=true;
		//----------------------------------------------------------------//

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el

		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) { IsPureInEL=true; }
		else { IsPureEL=true; }
		//--------------------------------------------------------------------------------------------------------------------------------//

		//Get true start/end point -----------------------------------------------------------------------//
		double true_endz=primary_truth_EndPosition_MC[2]; 
		double true_endy=primary_truth_EndPosition_MC[1]; 
		double true_endx=primary_truth_EndPosition_MC[0];

		double true_stz=primary_truth_StartPosition_MC[2];
		double true_sty=primary_truth_StartPosition_MC[1];
		double true_stx=primary_truth_StartPosition_MC[0];

		bool IsTrueEndOutside=false;
		if (true_endz<0.) IsTrueEndOutside=true;

		//Get reco info ----------------------------------------------------------------------------------//
		//Evt Classification -----------------------------------------------------------------------------//
		//signal -----------------------------//
		bool kinel=false;
		bool kel=false;
		//bool kmcs=false;
		if (IsBeamMatch) { //beam-match
			if (IsPureInEL) kinel=true;
			if (IsPureEL) kel=true;
			//if (IsPureMCS) kmcs=true;
		} //beam-match

		//background ------------------------------------------------------------------------//
		bool kMIDcosmic=false; //beam or cosmic
		bool kMIDpi=false; //+-pi
		bool kMIDp=false; //p
		bool kMIDmu=false; //mu
		bool kMIDeg=false; //e/gamma
		bool kMIDother=false; //other
		if (!IsBeamMatch) { //!beam-match
			if (primary_truth_byE_origin==2) { 
				kMIDcosmic=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==211) {
				kMIDpi=true;
			}
			else if (primary_truth_byE_PDG==2212) {
				kMIDp=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==13) {
				kMIDmu=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==11 || primary_truth_byE_PDG==22) {
				kMIDeg=true;
			}
			else {
				kMIDother=true;
			}
		} //!beam-match	
		//cout<<"kMIDcosmic:"<<kMIDcosmic<<endl;
		//Evt Classification -----------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		//double dx_reco_stx__bposx_ff=-999, dy_reco_sty__bposy_ff=-999, dz_reco_stz__bposz_ff=-999;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

			if (reco_stz>reco_endz) {
				reco_endx=primtrk_hitx->at(0); 
				reco_endy=primtrk_hity->at(0);
				reco_endz=primtrk_hitz->at(0);

				reco_stx=primtrk_hitx->at(primtrk_dedx->size()-1);	
				reco_sty=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_stz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}


			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	
			Fill1DHist(h1d_stxy_CaloSz, beam_dxy);

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx
		}

		//cosine_theta/cut
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999; 
		//cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction(no SCE corr.)
		TVector3 dir;
		if (IsCaloSize) {	
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			//TVector3 dir = pt1 - pt0;
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
      			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
      			TVector3 beamdir(beamDirx_spec->at(0),beamDiry_spec->at(0),beamDirz_spec->at(0));
      			beamdir = beamdir.Unit();
      			//beam_costh = dir.Dot(beamdir);
      			cosine_beam_spec_primtrk=dir.Dot(beamdir);
		}

		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

		//xy-cut
		bool IsXY=false;		
		double reco_stx_noSCE=0, reco_sty_noSCE=0, reco_stz_noSCE=0; //start-pos, before sce
		if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
			reco_stx_noSCE=primaryStartPosition[0];
			reco_sty_noSCE=primaryStartPosition[1];
			reco_stz_noSCE=primaryStartPosition[2];
		} //check if Pandora flip the sign
		else {
			reco_stx_noSCE=primaryEndPosition[0];
			reco_sty_noSCE=primaryEndPosition[1];
			reco_stz_noSCE=primaryEndPosition[2];
		}
		//if ((pow(((reco_stx_noSCE-mean_x)/dev_x),2)+pow(((reco_sty_noSCE-mean_y)/dev_y),2))<=1.) IsXY=true;

		//beam quality cut
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=-999;
		vector<double> reco_trklen_accum;
  		reco_trklen_accum.reserve(primtrk_hitz->size());
		double kereco_calo=0;
		double kereco_range=0;
		double kereco_range2=0;
		vector<double> EDept;
		double pid=-99;
		if (IsCaloSize) { //if calo size not empty
		  vector<double> trkdedx;
		  vector<double> trkres;
		  for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a given track
			double hitx_reco=primtrk_hitx->at(h);
			double hity_reco=primtrk_hity->at(h);
			double hitz_reco=primtrk_hitz->at(h);
			double resrange_reco=primtrk_resrange->at(h);

			double dqdx=primtrk_dqdx->at(h);
			double pitch=primtrk_pitch->at(h);

			int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
			double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

			if (wid_reco>wid_reco_max) { 
				wid_reco_max=wid_reco;
				index_reco_endz=(int)-1+primtrk_wid->size()-h;
			}

			double cali_dedx=0.;
			//cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
			//EDept.push_back(cali_dedx*pitch);

			if (h==1) range_reco=0;
			if (h>=1) {
    					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
					    		    pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
					    		    pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
			}

			kereco_calo+=cali_dedx*pitch;
			//kereco_range+=pitch*dedx_predict(resrange_reco);
			//kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

			trkdedx.push_back(cali_dedx);
			trkres.push_back(resrange_reco);

		  } //loop over reco hits of a given track
		  //range_reco=primtrk_range->at(0);

		  pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis
		} //if calo size not empty


		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		double bx_spec=beamPosx_spec->at(0);
		double by_spec=beamPosy_spec->at(0);

		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		//if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			//if (pid>pid_1) IsRecoInEL=true; 
			//if (pid<=pid_1) IsRecoStop=true; 
		//} //inel region
		//if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			//if (pid>pid_2) IsRecoInEL=true; 
			//if (pid<=pid_2) IsRecoStop=true;
		//} //stopping p region


		//kinetic energies
		//double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		//double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		//double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		//double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]
		//double ke_calo_MeV=0;

		//Get true trklen ---------------------------------------------------------------------------------------//
		double range_true=-999;
		int key_st = 0;
		double tmp_z = 9999;
		vector<double> true_trklen_accum;
		//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
		true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		//cout<<"ck0"<<endl;
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			//cout<<"ck0/"<<endl;
			true_trklen_accum[iz] = 0.; // initialize true_trklen_accum
			//cout<<"ck0///"<<endl;
		}
		//cout<<"ck1"<<endl;
		for (int iz=key_st+1; iz<(int)beamtrk_z->size(); iz++){
			if (iz == key_st+1) range_true = 0;
			range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
			true_trklen_accum[iz] = range_true;
		}




		if (IsPandoraSlice&&IsCaloSize) { //calosz cut
			//before bmrw
			Fill1DHist(h1d_trklen_CaloSz, range_reco);
			Fill1DHist(h1d_endz_CaloSz, reco_endz);
			Fill1DHist(h1d_stz_CaloSz, reco_stz);
			Fill1DHist(h1d_sty_CaloSz, reco_sty);
			Fill1DHist(h1d_stx_CaloSz, reco_stx);



/*
			Fill1DHist(h1d_truetrklen_CaloSz, range_true);
			Fill1DHist(h1d_trueendz_CaloSz, true_endz);

			if (kinel) Fill1DHist(h1d_trklen_CaloSz_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_CaloSz_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_CaloSz_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_CaloSz_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_CaloSz_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_CaloSz_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_CaloSz_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_CaloSz_midother, range_reco);
*/

		} //calosz cut


		if (IsPos&&IsPandoraSlice&&IsCaloSize) { //Pos
			Fill1DHist(h1d_trklen_Pos, range_reco);
			Fill1DHist(h1d_cosine_Pos, cosine_beam_spec_primtrk);
		} //Pos

		if (IsBQ&&IsPos&&IsPandoraSlice&&IsCaloSize) { //BQ
			Fill1DHist(h1d_trklen_BQ, range_reco);
			bx_by_RecoAll->Fill(bx_spec, by_spec);
			h1d_ntrklen_BQ->Fill(range_reco/csda_val_spec);
			h1d_pid_BQ->Fill(pid);
		} //BQ
	} //main entry loop

	//save results ---------------------------------------------------------//
   	//TFile *fout = new TFile("mc_sceoff.root","RECREATE");
   	TFile *fout = new TFile("mc_sceoff_new.root","RECREATE");

		h1d_trklen_CaloSz->Write();
		h1d_endz_CaloSz->Write();
		h1d_stz_CaloSz->Write();
		h1d_sty_CaloSz->Write();
		h1d_stx_CaloSz->Write();
		h1d_stxy_CaloSz->Write();
	
		h1d_cosine_Pos->Write();
		h1d_trklen_Pos->Write();

		h1d_trklen_BQ->Write();

		bx_by_RecoAll->Write();
 		h1d_ntrklen_BQ->Write();
		h1d_pid_BQ->Write();


/*
		h1d_trklen_CaloSz_inel->Write();
		h1d_trklen_CaloSz_el->Write();
		h1d_trklen_CaloSz_midcosmic->Write();
		h1d_trklen_CaloSz_midpi->Write();
		h1d_trklen_CaloSz_midp->Write();
		h1d_trklen_CaloSz_midmu->Write();
		h1d_trklen_CaloSz_mideg->Write();
		h1d_trklen_CaloSz_midother->Write();
*/

	fout->Close();



}