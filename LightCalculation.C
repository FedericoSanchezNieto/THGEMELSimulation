#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TRandom.h"
#include "TCanvas.h"
#include <iostream>

TRandom *r;

double pressure = 2.; // in bars 

double CalcEN(double field = 0, double p = 0){
    // gas parameteres
    double n(0.);
    double N0(6.022e23);
    double R(8.3145);//J/mol/K
    double T(273.);//K
    //Geometry
    double cm_unit = 1e-6;
    double toTd = 1e-17;
    n = p*N0/(R*T)*cm_unit;
    // std::cout << n << std::endl;
    return field/n/toTd;
}


double photons(double E, double d) {
  double p = pressure*100000.; // Pascal 
 
  double T = CalcEN(E/100.,p);

  double T0 = 5.0; // 4.5;
  double Yield = 18.8;

  d *=100.;
 
  if( T > T0 ) 
    return Yield*(T - T0)*d ;
  else
    return 0; 
}


double amplification(double E, double d) {

  if( CalcEN(E,pressure*100000.) < 12. ) return 1.; // Ionisation threshold
  
  double p = pressure*750.;  // Torr

  E /= 100.;  // in V/cm 
  
  d *= 100.;  // d in cm 

  // Numbers from CERN PPE/93-50 
  
  double A = 2.8 * p * d ;  
  
  double B = 60.6; 
  
  double alpha = A*TMath::Exp(- B * p/E);

  // std::cout << E/p << " alpha " << alpha << std::endl; 
  
  //    std::cout << " E field " << E << "V/cm  new electrons: " << alpha << std::endl; 


  //  alpha *= 0.25; 
  
  return TMath::Exp(alpha); 					      
}


bool checkendinthgem(double x,double y) {

  double THGEMxmin  = 0.0;
  double THGEMxmax  = 0.0008; 

  double holeymin = 0.0490;
  double holeymax = 0.0494;

  if( x <= THGEMxmax && x >= THGEMxmin ) 
    if( y <= holeymin || y >=holeymax ) return true;

  return false;
}

#define NTRIESPEN 100 

double PENacceptance(double x,double dx,double dy,double dz){
    
  double hexradious = 5.6e-3; // in m   
  double hexwidth = 0.5e-3;
  double hexthick = 0.75e-3 ; 

  double rmin = hexradious; 
  double rmax = hexradious+hexthick; 

  
  // Assume intersection with a plane at x = 0; 
  // Then, smear the position radially and in phi in the plane, assume circles and not hexagones.  
 
  double nphotons = 0.0;

  double dd = hexwidth/dx; // hexwidth+x = x + dd * cos 
  
  for(int n = 0; n <  NTRIESPEN; n++ ) {
    double rad = TMath::Sqrt(r->Uniform(0.,rmax*rmax));
    if( rad > rmin ) continue;  // Photon hit the upper part of the mesh 
   
    double phi = r->Uniform(0.,2.*3.141592); 
    double y0 = rad*TMath::Cos(phi); 
    double z0 = rad*TMath::Sin(phi);  

    double y = y0+dd*dy;
    double z = z0+dd*dz; 

    if( y*y+z*z > rmin*rmin ) continue; // Hit bewteen the upper and the lower part of the mesh 

    nphotons++;  // Photon survived. 
    
  }
  
  return nphotons/(double)NTRIESPEN;

#if 0

  double THGEMPENPosition = 0.0034;
  double PENLength = 0.1;
  double THGEMlength= 0.1;
  
  double nphotons = 0;
  
  for(int n = 0; n <  NTRIESPEN; n++ ) {
    double x0 = x; 
    double y0 = r->Uniform(-THGEMlength/2.,THGEMlength/2.); 
    double z0 = r->Uniform(-THGEMlength/2.,THGEMlength/2.);

    double dd = (THGEMPENPosition-x0)/dx;

    double y = y0+dd*dy;
    if( TMath::Abs(y) > PENLength/2.) continue; 
    
    double z = z0+dd*dz; 
    if( TMath::Abs(z) > PENLength/2.) continue; 

    nphotons++; 
    
  }
  
  return nphotons/(double)NTRIESPEN;
#endif 
}


#define NTRIES 1000 

double acceptance(double x, double y) {
  
  double thgemxmin  = 0.0;
  double thgemxmax  = 0.0008; 

  double nphotons = 0;
  
  if( x > thgemxmax ) {

    for(int i = 0; i < NTRIES; i++) {
      double cos = r->Uniform(0.,1.);
      double phi = r->Uniform(-3.141592,3.141592);
      double sin = TMath::Sqrt(1.-cos*cos); 
      double dx = cos;
      double dy = sin*TMath::Cos(phi);
      double dz = sin*TMath::Sin(phi);
      nphotons+= PENacceptance(x,dx,dy,dz);
    }
  }
  else { 
    double holeymin = 0.0490;
    double holeymax = 0.0494;
    
    double ycenter = (holeymin+holeymax)/2.;
    double rad = (holeymax-holeymin)/2.; 
    double r2 = rad*rad; 
    
    // We consider a 3D but from the 2D. We build it by rotating the system and aligning the electron to the radious at z = 0 
    double zcenter = 0; 
    
    for(int i = 0; i < NTRIES; i++) {
      double cos = r->Uniform(0.,1.);
      double phi = r->Uniform(-3.141592,3.141592);
      double sin = TMath::Sqrt(1.-cos*cos); 
      double dx = cos;
      double dy = sin*TMath::Cos(phi);
      double dz = sin*TMath::Sin(phi);

      double z = zcenter;
      
      //First entry hole for photons before the THGEM 
      if( x < thgemxmin ) {
	double ddmin = (thgemxmin-x)/dx;  // Should be always positive 
	double yp1 = y + ddmin * dy;
	double zp1 = z + ddmin * dz;  
	if( (yp1-ycenter)*(yp1-ycenter)+(zp1-zcenter)*(zp1-zcenter) > r2 ) continue;
	//      std::cout << (yp1-ycenter)*(yp1-ycenter)+(zp1-zcenter)*(zp1-zcenter) << "  " <<  r2 << std::endl;
      }
      
      double ddmax = (thgemxmax-x)/dx;  // Should be always positive 
      double yp2 = y + ddmax * dy;
      double zp2 = z + ddmax * dz;  
      if( (yp2-ycenter)*(yp2-ycenter)+(zp2-zcenter)*(zp2-zcenter) > r2 ) continue;
      
      nphotons+= PENacceptance(x,dx,dy,dz); 
    
    }
  }
  
  double acceptVal  = nphotons/(double)NTRIES;  
  
  return 0.5*acceptVal; // Only 50% of the photons go forward. 
}


double LightCalculation(string file, int nelectrons = 1, double delta = 1.e-5, string value = "" , bool debug = false) {

  TFile *_file0 = TFile::Open(file.c_str());

  if( !_file0 ) { std::cout << " File " << file << " not found " << std::endl; return 0.;} 

  TH2D *Ex = (TH2D*) gROOT->FindObject("ex");

  TH2D *Ey = (TH2D*) gROOT->FindObject("ey"); 


  TObject *obj; 
  
  if( (obj = gROOT->FindObject("fieldplots")) ) delete obj;
  TCanvas *fieldPlots = new TCanvas("fieldplots","",800,800);
  fieldPlots->Divide(1,2);
  fieldPlots->cd(1); 
  Ex->Draw("colz");
  fieldPlots->cd(2); 
  Ey->Draw("colz");
  fieldPlots->Update(); 
  
  double xmin = Ex->GetXaxis()->GetXmin();
  double xmax = Ex->GetXaxis()->GetXmax(); 

  double ymin = Ex->GetYaxis()->GetXmin(); 
  double ymax = Ex->GetYaxis()->GetXmax(); 

  int maxit = 1e+5;

  if( debug ) {
    std::cout << " X = ( " << xmin << "," << xmax << " ) " << std::endl;
    std::cout << " Y = ( " << ymin << "," << ymax << " ) " << std::endl;
  }

  if( (obj = gROOT->FindObject("Photons")) ) delete obj;
  if( (obj = gROOT->FindObject("PhotonsAccp")) ) delete obj;
  if( (obj = gROOT->FindObject("PhotonsNoAmplAccp")) ) delete obj;
  if( (obj = gROOT->FindObject("Electrons")) ) delete obj;
  
  TH2D *Photons = new TH2D("Photons"," ",500,xmin,xmax,500,ymin,ymax);
  TH2D *PhotonsAccp = new TH2D("PhotonsAccp"," ",500,xmin,xmax,500,ymin,ymax);
  TH2D *PhotonsNoAmplAccp = new TH2D("PhotonsNoAmplAccp"," ",500,xmin,xmax,500,ymin,ymax);
  
  TH2D *Electrons = new TH2D("Electrons"," ",500,xmin,xmax,500,ymin,ymax);

  
  r = new TRandom;
  
  // Transport Electrons

  double totelectrons = 0;
  
  for(int i = 0; i < nelectrons; i++ ) {
    
    // Initiate the coordinates

#if 0     
    double x = xmin + 2.e-4;
    double y = r->Uniform(ymin,ymax); 
#else
    double x = -1.e-5;
    double y = r->Uniform(0.0494,ymax); 
#endif 
    
    bool alive = true; 

    int it = 0; 

#define MAXINXTRJ 10000
    
    double trajx[MAXINXTRJ];
    double trajy[MAXINXTRJ];

    int indxtrj = 0; 

    double electrons = 1.;


    double eprev = 0; 
    
    while(alive){
    
      double ex = Ex->Interpolate(x,y);
      double ey = Ey->Interpolate(x,y); 
     
      double e = TMath::Sqrt(ex*ex+ey*ey);
      
      double dx = -ex/e;
      double dy = -ey/e;

      if( eprev != 0. ) {double e1 = e; e = (e+eprev)/2.; eprev = e1; } //Take average in the transport  
      
      if( debug ) 
	std::cout << " Coord = ( " << x << ", " << y <<  "),  E = ( " <<  ex << " , " << ey << ") , d = (  " << dx << " ,  " << dy << " ) " << std::endl; 
      
      // new coordinate:
      
      x += dx*delta;
      y += dy*delta;

      if( x > xmax || x < xmin || y > ymax || y < ymin || it > maxit ) {alive= false; break;}

      if( checkendinthgem(x,y) ) { alive = false;break;}

      double ampl = amplification(e,delta); // std::cout << " Electrons " << electrons << std::endl;  

      Electrons->Fill(x,y,ampl-1.); 
      
      electrons *= ampl;

      double nphotons = photons(e,delta)*electrons;
      double nphotonsNoAmpl = photons(e,delta); 

      if( nphotons > 0 ) {
      
	Photons->Fill(x,y,nphotons);

	double acc = acceptance(x,y);

	PhotonsAccp->Fill(x,y,acc*nphotons); 

	PhotonsNoAmplAccp->Fill(x,y,nphotonsNoAmpl*acc); 

      } 

	
      if( indxtrj < MAXINXTRJ ) {
	trajx[indxtrj] = x;
	trajy[indxtrj] = y; 
	indxtrj++;
      }
      
      it++; 
    }

    totelectrons += electrons; 
    
    if( i < 100 ) {
    
      TGraph *electron = new TGraph(indxtrj,trajx,trajy); 

      fieldPlots->cd(1); 
      electron->Draw("same");
      fieldPlots->Update();
      
      fieldPlots->cd(2); 
      electron->Draw("same");
      fieldPlots->Update();
    }
  }

  TH2D *Accp = (TH2D*) PhotonsAccp->Clone();
  Accp->Divide(Photons); 


  if( (obj = gROOT->FindObject("Photonplots")) ) delete obj;
  TCanvas *PhotonPlots = new TCanvas("Photonplots","",800,800);
  PhotonPlots->Divide(1,3);
  PhotonPlots->cd(1); 
  Photons->Draw("colz");
  PhotonPlots->cd(2);
  PhotonsAccp->Draw("colz");
  PhotonPlots->cd(3);
  Accp->Draw("colz"); 
  PhotonPlots->Update();


  if( (obj = gROOT->FindObject("Electronplots")) ) delete obj;
  TCanvas *ElectronPlots = new TCanvas("Electronplots","",800,800);
  Electrons->Draw("colz");
  ElectronPlots->Update();

  
  std::cout<<  value << "  " <<  PhotonsNoAmplAccp->Integral()/(double)nelectrons  << "   "  << PhotonsAccp->Integral()/(double)nelectrons << "  "  << totelectrons/(double) nelectrons  <<  std::endl;

  // _file0->Close(); 

  return   PhotonsAccp->Integral()/(double)nelectrons; 
}


void RunAll(int nel = 10, double P = 2. ){

  pressure = P; 
  
  LightCalculation("Histos_meshNeg400V.root",nel,1.e-5,"-400");
  LightCalculation("Histos_meshNeg300V.root",nel,1.e-5,"-300");
  LightCalculation("Histos_meshNeg200V.root",nel,1.e-5,"-200");
  LightCalculation("Histos_meshNeg100V.root",nel,1.e-5,"-100");
  LightCalculation("Histos_mesh0V.root",nel,1.e-5,"0");
  LightCalculation("Histos_mesh100V.root",nel,1.e-5,"100");
  LightCalculation("Histos_mesh150V.root",nel,1.e-5,"150");
  LightCalculation("Histos_mesh200V.root",nel,1.e-5,"200");
  LightCalculation("Histos_mesh300V.root",nel,1.e-5,"300");
  LightCalculation("Histos_mesh500V.root",nel,1.e-5,"500");
  LightCalculation("Histos_mesh600V.root",nel,1.e-5,"600");
  LightCalculation("Histos_mesh900V.root",nel,1.e-5,"900");
}
