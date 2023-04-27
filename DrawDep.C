void DrawDep(string filename) {
  std::ifstream myfile;
  myfile.open(filename);
  std::string myline;
  double varr[20];
  double nogainarr[20];
  double gainarr[20];
  double elarr[20];
  int indx = 0;
  double max = 0;
  double vmin = 0;
  double vmax = 0; 
  
  if ( myfile.is_open() ) {
    while ( myfile ) { // equivalent to myfile.good()
      double v,gain,nogain,el; 
      std::getline (myfile, myline);
      sscanf(myline.c_str()," %lg %lg %lg %lg ",&v,&nogain,&gain,&el);
      varr[indx] = v;  
      nogainarr[indx] = nogain;
      gainarr[indx] = gain;
      elarr[indx] = el;
      if( el > max ) max = el;
      if( nogain > max ) max = nogain;
      if( gain > max ) max = gain;

      if( vmin > v ) vmin =v;
      if( vmax < v ) vmax =v;
      
      indx++;
    }
    TGraph *ggain = new TGraph(indx,varr,gainarr);
    TGraph *gnogain = new TGraph(indx,varr,nogainarr);
    TGraph *gel = new TGraph(indx,varr,elarr); 

    TH2F *h2 = new TH2F("h","",100,vmin-10,vmax+10,100,0,max); 
    h2->Draw(); 
    ggain->Draw("LP");
    ggain->SetLineColor(2); 
    gnogain->Draw("LP same");
    gnogain->SetLineColor(3); 
    gel->Draw("LP same");
  }
  else {
    std::cout << "Couldn't open file\n";
  }
  
}
