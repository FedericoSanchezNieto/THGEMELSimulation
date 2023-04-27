void CreateNtuple(string filein, string fileout) {

  TNtuple * n = new TNtuple("Field","Field","x:y:E:Ex:Ey:V");
  
  n->ReadFile(filein.c_str());
  
  n->SaveAs(fileout.c_str()); 
   
}
 
