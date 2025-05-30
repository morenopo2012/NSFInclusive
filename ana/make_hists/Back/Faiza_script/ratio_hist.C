void ratio_hist(){


TH1D* h1 = new TH1D("h1","h1",13,0.,13.); //Iron: MENNDL
   h1->SetBinContent(7,0.7380714);
   h1->SetBinContent(8,0.3461664);
   h1->SetBinContent(9,3.307287);
   h1->SetBinContent(10,14.75906);
   h1->SetBinContent(11,447.1312);
   h1->SetBinContent(12,8.649494);
   h1->SetBinContent(13,2.0808);
   h1->SetBinError(7,0.5225792);
   h1->SetBinError(8,0.3461664);
   h1->SetBinError(9,1.102973);
   h1->SetBinError(10,2.312106);
   h1->SetBinError(11,12.76682);
   h1->SetBinError(12,1.767272);
   h1->SetBinError(13,0.850575);

TH1D* h2 = new TH1D("h2","h2",13,0.,13.); //Iron: Caffe
   h2->SetBinContent(7,0.3501377);
   h2->SetBinContent(8,0.3511236);
   h2->SetBinContent(9,5.398327);
   h2->SetBinContent(10,12.75721);
   h2->SetBinContent(11,458.6205);
   h2->SetBinContent(12,3.956676);
   h2->SetBinContent(13,1.742695);
   h2->SetBinError(7,0.3501377);
   h2->SetBinError(8,0.3511236);
   h2->SetBinError(9,1.394452);
   h2->SetBinError(10,2.161792);
   h2->SetBinError(11,12.92806);
   h2->SetBinError(12,1.194632);
   h2->SetBinError(13,0.7794458);
 
   TH1D* hratio = (TH1D*)h1->Clone("hratio");
   hratio->Divide(h1,h2);

   for(int i = 0; i < hratio->GetNbinsX(); i++){
      cout<<"h2->SetBinContent("<<i+1<<","<<hratio->GetBinContent(i+1)<<");"<<endl;
      cout<<"h2->SetBinError("  <<i+1<<","<<hratio->GetBinError(i+1)<<");"<<endl;
   }

   hratio->Draw("ep");

}
