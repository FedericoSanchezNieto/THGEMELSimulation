void MakeHistos(string name){
  double xmin = -0.0010;
  double xmax = 0.0034; 

  double ymin = 0.0488; 
  double ymax = 0.0496; 

  
  TH2D *ex = new TH2D("ex"," ",500,xmin,xmax,500,ymin,ymax); 
  TH2D *ey = new TH2D("ey"," ",500,xmin,xmax,500,ymin,ymax); 
  TH2D *n = new TH2D("n"," ",500,xmin,xmax,500,ymin,ymax); 

  Field->Project("ex","y:x","Ex","colz");
  Field->Project("ey","y:x","Ey","colz");
  Field->Project("n","y:x","","colz");

  ex->Divide(n);
  ey->Divide(n);

  ex->Draw("colz");

  TFile *f = new TFile(name.c_str(),"RECREATE");

  ex->Write();
  ey->Write();

  f->Close(); 
  
}
