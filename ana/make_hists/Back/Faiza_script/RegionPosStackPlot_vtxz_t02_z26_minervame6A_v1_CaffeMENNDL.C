{
//=========Macro generated from canvas: RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL/RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL
//=========  (Wed Oct 30 23:53:52 2019) by ROOT version5.34/36
   TCanvas *RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL = new TCanvas("RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL", "RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL",0,0,1200,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetHighLightColor(2);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->Range(448.1211,-118.4211,492.336,671.0526);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetFillColor(0);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetBorderMode(0);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetBorderSize(2);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetGridx();
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetGridy();
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetLeftMargin(0.15);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetRightMargin(0.15);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetBottomMargin(0.15);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetFrameLineWidth(2);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetFrameBorderMode(0);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetFrameLineWidth(2);
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetFrameBorderMode(0);
   Double_t xAxis16[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0000_4962 = new TH1D("tmpMC_0000_4962","Carbon",92, xAxis16);
   tmpMC_0000_4962->SetMinimum(0);
   tmpMC_0000_4962->SetMaximum(600);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetFillColor(ci);
   tmpMC_0000_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetLineColor(ci);
   tmpMC_0000_4962->SetLineWidth(5);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetMarkerColor(ci);
   tmpMC_0000_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0000_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0000_4962->GetXaxis()->CenterTitle(true);
   tmpMC_0000_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetXaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetXaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0000_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetYaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetYaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleFont(42);
   tmpMC_0000_4962->Draw("HIST");
   
   THStack *hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("Stacked 1D histograms");
   Double_t xAxis17[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1F *hs_stack_3 = new TH1F("hs_stack_3","Stacked 1D histograms",92, xAxis17);
   hs_stack_3->SetMinimum(0);
   hs_stack_3->SetMaximum(598.1129);
   hs_stack_3->SetDirectory(0);
   hs_stack_3->SetStats(0);

   ci = TColor::GetColor("#000099");
   hs_stack_3->SetLineColor(ci);
   hs_stack_3->SetLineWidth(2);
   hs_stack_3->SetMarkerStyle(20);
   hs_stack_3->SetMarkerSize(1.75);
   hs_stack_3->GetXaxis()->SetLabelFont(42);
   hs_stack_3->GetXaxis()->SetLabelSize(0.05);
   hs_stack_3->GetXaxis()->SetTitleSize(0.06);
   hs_stack_3->GetXaxis()->SetTitleOffset(1.15);
   hs_stack_3->GetYaxis()->SetLabelFont(42);
   hs_stack_3->GetYaxis()->SetLabelSize(0.05);
   hs_stack_3->GetYaxis()->SetTitleSize(0.06);
   hs_stack_3->GetYaxis()->SetTitleOffset(1.2);
   hs_stack_3->GetZaxis()->SetLabelFont(42);
   hs_stack_3->GetZaxis()->SetLabelSize(0.05);
   hs_stack_3->GetZaxis()->SetTitleSize(0.06);
   hs_stack_3->GetZaxis()->SetTitleOffset(0.75);
   hs->SetHistogram(hs_stack_3);
   
   Double_t xAxis18[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0000_4962 = new TH1D("tmpMC_0000_4962","Carbon",92, xAxis18);
   tmpMC_0000_4962->SetMinimum(0);
   tmpMC_0000_4962->SetMaximum(600);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetFillColor(ci);
   tmpMC_0000_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetLineColor(ci);
   tmpMC_0000_4962->SetLineWidth(5);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetMarkerColor(ci);
   tmpMC_0000_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0000_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0000_4962->GetXaxis()->CenterTitle(true);
   tmpMC_0000_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetXaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetXaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0000_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetYaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetYaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0000_4962,"");
   Double_t xAxis19[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0001_4962 = new TH1D("tmpMC_0001_4962","Lead",92, xAxis19);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetFillColor(ci);
   tmpMC_0001_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetLineColor(ci);
   tmpMC_0001_4962->SetLineWidth(5);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetMarkerColor(ci);
   tmpMC_0001_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0001_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0001_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0001_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0001_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0001_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0001_4962,"");
   Double_t xAxis20[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0002_4962 = new TH1D("tmpMC_0002_4962","Iron",92, xAxis20);
   tmpMC_0002_4962->SetBinContent(7,0.7380714);
   tmpMC_0002_4962->SetBinContent(8,0.3461664);
   tmpMC_0002_4962->SetBinContent(9,3.307287);
   tmpMC_0002_4962->SetBinContent(10,14.75906);
   tmpMC_0002_4962->SetBinContent(11,447.1312);
   tmpMC_0002_4962->SetBinContent(12,8.649494);
   tmpMC_0002_4962->SetBinContent(13,2.0808);
   tmpMC_0002_4962->SetBinError(7,0.5225792);
   tmpMC_0002_4962->SetBinError(8,0.3461664);
   tmpMC_0002_4962->SetBinError(9,1.102973);
   tmpMC_0002_4962->SetBinError(10,2.312106);
   tmpMC_0002_4962->SetBinError(11,12.76682);
   tmpMC_0002_4962->SetBinError(12,1.767272);
   tmpMC_0002_4962->SetBinError(13,0.850575);
   tmpMC_0002_4962->SetEntries(1319);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetFillColor(ci);
   tmpMC_0002_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetLineColor(ci);
   tmpMC_0002_4962->SetLineWidth(5);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetMarkerColor(ci);
   tmpMC_0002_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0002_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0002_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0002_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0002_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0002_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0002_4962,"");
   Double_t xAxis21[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0003_4962 = new TH1D("tmpMC_0003_4962","Upstream Scint.",92, xAxis21);
   tmpMC_0003_4962->SetBinContent(7,54.17861);
   tmpMC_0003_4962->SetBinContent(8,94.64936);
   tmpMC_0003_4962->SetBinContent(9,91.53072);
   tmpMC_0003_4962->SetBinContent(10,68.08932);
   tmpMC_0003_4962->SetBinContent(11,44.08771);
   tmpMC_0003_4962->SetBinContent(12,0.3817845);
   tmpMC_0003_4962->SetBinContent(14,0.3837732);
   tmpMC_0003_4962->SetBinError(7,4.496386);
   tmpMC_0003_4962->SetBinError(8,5.87842);
   tmpMC_0003_4962->SetBinError(9,5.780719);
   tmpMC_0003_4962->SetBinError(10,4.932879);
   tmpMC_0003_4962->SetBinError(11,3.997932);
   tmpMC_0003_4962->SetBinError(12,0.3817845);
   tmpMC_0003_4962->SetBinError(14,0.3837732);
   tmpMC_0003_4962->SetEntries(977);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetFillColor(ci);
   tmpMC_0003_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetLineColor(ci);
   tmpMC_0003_4962->SetLineWidth(5);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetMarkerColor(ci);
   tmpMC_0003_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0003_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0003_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0003_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0003_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0003_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0003_4962,"");
   Double_t xAxis22[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0004_4962 = new TH1D("tmpMC_0004_4962","Downstream Scint.",92, xAxis22);
   tmpMC_0004_4962->SetBinContent(7,0.3437299);
   tmpMC_0004_4962->SetBinContent(8,1.40444);
   tmpMC_0004_4962->SetBinContent(9,0.3437158);
   tmpMC_0004_4962->SetBinContent(10,0.3599438);
   tmpMC_0004_4962->SetBinContent(11,65.79146);
   tmpMC_0004_4962->SetBinContent(12,48.4214);
   tmpMC_0004_4962->SetBinContent(13,110.4223);
   tmpMC_0004_4962->SetBinContent(14,100.9683);
   tmpMC_0004_4962->SetBinContent(15,72.76481);
   tmpMC_0004_4962->SetBinError(7,0.3437299);
   tmpMC_0004_4962->SetBinError(8,0.7022219);
   tmpMC_0004_4962->SetBinError(9,0.3437158);
   tmpMC_0004_4962->SetBinError(10,0.3599438);
   tmpMC_0004_4962->SetBinError(11,4.96303);
   tmpMC_0004_4962->SetBinError(12,4.187757);
   tmpMC_0004_4962->SetBinError(13,6.341359);
   tmpMC_0004_4962->SetBinError(14,6.043206);
   tmpMC_0004_4962->SetBinError(15,5.143353);
   tmpMC_0004_4962->SetEntries(1108);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetFillColor(ci);
   tmpMC_0004_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetLineColor(ci);
   tmpMC_0004_4962->SetLineWidth(5);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetMarkerColor(ci);
   tmpMC_0004_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0004_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0004_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0004_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0004_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0004_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0004_4962,"");
   Double_t xAxis23[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0005_4962 = new TH1D("tmpMC_0005_4962","Other Target",92, xAxis23);
   tmpMC_0005_4962->SetBinContent(7,5.412083);
   tmpMC_0005_4962->SetBinContent(8,2.854513);
   tmpMC_0005_4962->SetBinContent(9,3.988622);
   tmpMC_0005_4962->SetBinContent(10,0.7407119);
   tmpMC_0005_4962->SetBinContent(11,12.62097);
   tmpMC_0005_4962->SetBinContent(12,2.831497);
   tmpMC_0005_4962->SetBinContent(13,4.546255);
   tmpMC_0005_4962->SetBinContent(14,9.24507);
   tmpMC_0005_4962->SetBinContent(15,17.18615);
   tmpMC_0005_4962->SetBinError(7,1.398722);
   tmpMC_0005_4962->SetBinError(8,1.010255);
   tmpMC_0005_4962->SetBinError(9,1.204079);
   tmpMC_0005_4962->SetBinError(10,0.524135);
   tmpMC_0005_4962->SetBinError(11,2.135405);
   tmpMC_0005_4962->SetBinError(12,1.001729);
   tmpMC_0005_4962->SetBinError(13,1.322914);
   tmpMC_0005_4962->SetBinError(14,1.856676);
   tmpMC_0005_4962->SetBinError(15,2.482675);
   tmpMC_0005_4962->SetEntries(164);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetFillColor(ci);
   tmpMC_0005_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetLineColor(ci);
   tmpMC_0005_4962->SetLineWidth(5);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetMarkerColor(ci);
   tmpMC_0005_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0005_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0005_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0005_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0005_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0005_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0005_4962,"");
   hs->Draw("same hist");
   
   hs = new THStack();
   hs->SetName("hs");
   hs->SetTitle("Stacked 1D histograms");
   Double_t xAxis24[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1F *hs_stack_3_stack_4 = new TH1F("hs_stack_3_stack_4","Stacked 1D histograms",92, xAxis24);
   hs_stack_3_stack_4->SetMinimum(0);
   hs_stack_3_stack_4->SetMaximum(598.1129);
   hs_stack_3_stack_4->SetDirectory(0);
   hs_stack_3_stack_4->SetStats(0);

   ci = TColor::GetColor("#000099");
   hs_stack_3_stack_4->SetLineColor(ci);
   hs_stack_3_stack_4->SetLineWidth(2);
   hs_stack_3_stack_4->SetMarkerStyle(20);
   hs_stack_3_stack_4->SetMarkerSize(1.75);
   hs_stack_3_stack_4->GetXaxis()->SetLabelFont(42);
   hs_stack_3_stack_4->GetXaxis()->SetLabelSize(0.05);
   hs_stack_3_stack_4->GetXaxis()->SetTitleSize(0.06);
   hs_stack_3_stack_4->GetXaxis()->SetTitleOffset(1.15);
   hs_stack_3_stack_4->GetYaxis()->SetLabelFont(42);
   hs_stack_3_stack_4->GetYaxis()->SetLabelSize(0.05);
   hs_stack_3_stack_4->GetYaxis()->SetTitleSize(0.06);
   hs_stack_3_stack_4->GetYaxis()->SetTitleOffset(1.2);
   hs_stack_3_stack_4->GetZaxis()->SetLabelFont(42);
   hs_stack_3_stack_4->GetZaxis()->SetLabelSize(0.05);
   hs_stack_3_stack_4->GetZaxis()->SetTitleSize(0.06);
   hs_stack_3_stack_4->GetZaxis()->SetTitleOffset(0.75);
   hs->SetHistogram(hs_stack_3_stack_4);
   
   Double_t xAxis25[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0000_4962 = new TH1D("tmpMC_0000_4962","Carbon",92, xAxis25);
   tmpMC_0000_4962->SetMinimum(0);
   tmpMC_0000_4962->SetMaximum(600);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetFillColor(ci);
   tmpMC_0000_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetLineColor(ci);
   tmpMC_0000_4962->SetLineWidth(5);

   ci = TColor::GetColor("#660066");
   tmpMC_0000_4962->SetMarkerColor(ci);
   tmpMC_0000_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0000_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0000_4962->GetXaxis()->CenterTitle(true);
   tmpMC_0000_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetXaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetXaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0000_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetYaxis()->SetLabelSize(0.05);
   tmpMC_0000_4962->GetYaxis()->SetTitleSize(0.06);
   tmpMC_0000_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0000_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0000_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0000_4962,"");
   Double_t xAxis26[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0001_4962 = new TH1D("tmpMC_0001_4962","Lead",92, xAxis26);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetFillColor(ci);
   tmpMC_0001_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetLineColor(ci);
   tmpMC_0001_4962->SetLineWidth(5);

   ci = TColor::GetColor("#006699");
   tmpMC_0001_4962->SetMarkerColor(ci);
   tmpMC_0001_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0001_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0001_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0001_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0001_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0001_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0001_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0001_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0001_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0001_4962,"");
   Double_t xAxis27[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0002_4962 = new TH1D("tmpMC_0002_4962","Iron",92, xAxis27);
   tmpMC_0002_4962->SetBinContent(7,0.7380714);
   tmpMC_0002_4962->SetBinContent(8,0.3461664);
   tmpMC_0002_4962->SetBinContent(9,3.307287);
   tmpMC_0002_4962->SetBinContent(10,14.75906);
   tmpMC_0002_4962->SetBinContent(11,447.1312);
   tmpMC_0002_4962->SetBinContent(12,8.649494);
   tmpMC_0002_4962->SetBinContent(13,2.0808);
   tmpMC_0002_4962->SetBinError(7,0.5225792);
   tmpMC_0002_4962->SetBinError(8,0.3461664);
   tmpMC_0002_4962->SetBinError(9,1.102973);
   tmpMC_0002_4962->SetBinError(10,2.312106);
   tmpMC_0002_4962->SetBinError(11,12.76682);
   tmpMC_0002_4962->SetBinError(12,1.767272);
   tmpMC_0002_4962->SetBinError(13,0.850575);
   tmpMC_0002_4962->SetEntries(1319);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetFillColor(ci);
   tmpMC_0002_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetLineColor(ci);
   tmpMC_0002_4962->SetLineWidth(5);

   ci = TColor::GetColor("#cc0000");
   tmpMC_0002_4962->SetMarkerColor(ci);
   tmpMC_0002_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0002_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0002_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0002_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0002_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0002_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0002_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0002_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0002_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0002_4962,"");
   Double_t xAxis28[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0003_4962 = new TH1D("tmpMC_0003_4962","Upstream Scint.",92, xAxis28);
   tmpMC_0003_4962->SetBinContent(7,54.17861);
   tmpMC_0003_4962->SetBinContent(8,94.64936);
   tmpMC_0003_4962->SetBinContent(9,91.53072);
   tmpMC_0003_4962->SetBinContent(10,68.08932);
   tmpMC_0003_4962->SetBinContent(11,44.08771);
   tmpMC_0003_4962->SetBinContent(12,0.3817845);
   tmpMC_0003_4962->SetBinContent(14,0.3837732);
   tmpMC_0003_4962->SetBinError(7,4.496386);
   tmpMC_0003_4962->SetBinError(8,5.87842);
   tmpMC_0003_4962->SetBinError(9,5.780719);
   tmpMC_0003_4962->SetBinError(10,4.932879);
   tmpMC_0003_4962->SetBinError(11,3.997932);
   tmpMC_0003_4962->SetBinError(12,0.3817845);
   tmpMC_0003_4962->SetBinError(14,0.3837732);
   tmpMC_0003_4962->SetEntries(977);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetFillColor(ci);
   tmpMC_0003_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetLineColor(ci);
   tmpMC_0003_4962->SetLineWidth(5);

   ci = TColor::GetColor("#cc00cc");
   tmpMC_0003_4962->SetMarkerColor(ci);
   tmpMC_0003_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0003_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0003_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0003_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0003_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0003_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0003_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0003_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0003_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0003_4962,"");
   Double_t xAxis29[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0004_4962 = new TH1D("tmpMC_0004_4962","Downstream Scint.",92, xAxis29);
   tmpMC_0004_4962->SetBinContent(7,0.3437299);
   tmpMC_0004_4962->SetBinContent(8,1.40444);
   tmpMC_0004_4962->SetBinContent(9,0.3437158);
   tmpMC_0004_4962->SetBinContent(10,0.3599438);
   tmpMC_0004_4962->SetBinContent(11,65.79146);
   tmpMC_0004_4962->SetBinContent(12,48.4214);
   tmpMC_0004_4962->SetBinContent(13,110.4223);
   tmpMC_0004_4962->SetBinContent(14,100.9683);
   tmpMC_0004_4962->SetBinContent(15,72.76481);
   tmpMC_0004_4962->SetBinError(7,0.3437299);
   tmpMC_0004_4962->SetBinError(8,0.7022219);
   tmpMC_0004_4962->SetBinError(9,0.3437158);
   tmpMC_0004_4962->SetBinError(10,0.3599438);
   tmpMC_0004_4962->SetBinError(11,4.96303);
   tmpMC_0004_4962->SetBinError(12,4.187757);
   tmpMC_0004_4962->SetBinError(13,6.341359);
   tmpMC_0004_4962->SetBinError(14,6.043206);
   tmpMC_0004_4962->SetBinError(15,5.143353);
   tmpMC_0004_4962->SetEntries(1108);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetFillColor(ci);
   tmpMC_0004_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetLineColor(ci);
   tmpMC_0004_4962->SetLineWidth(5);

   ci = TColor::GetColor("#ff9933");
   tmpMC_0004_4962->SetMarkerColor(ci);
   tmpMC_0004_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0004_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0004_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0004_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0004_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0004_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0004_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0004_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0004_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0004_4962,"");
   Double_t xAxis30[93] = {424.1122, 428.2245, 432.646, 437.0674, 441.4889, 445.9104, 450.3319, 454.7534, 459.1748, 463.5963, 468.0178, 472.4393, 476.8608, 481.2822, 485.7037, 490.1252, 494.5467, 498.9682, 503.3897, 507.8111, 512.2326, 516.4434, 544.5941, 549.0156, 553.4371, 557.8585, 562.28, 566.7015, 571.123, 575.5445, 579.9152, 584.4381, 588.961, 593.4839, 598.0068, 602.5297, 607.0526, 611.5755, 616.0984, 620.6213, 625.1442, 629.6671, 634.19, 638.7129, 643.2358, 647.7587, 652.2816, 656.8045, 661.3274, 665.8503, 670.3732, 674.8961, 679.419, 683.9419, 688.4648, 692.9877, 697.5106, 702.0335, 706.5564, 711.0793, 715.6022, 720.1251, 724.648, 729.1709, 733.6938, 738.2167, 742.7396, 747.2625, 751.7854, 756.3083, 760.8312, 765.3541, 769.877, 774.3999, 778.9228, 783.4457, 787.9686, 792.4915, 797.0144, 801.5373, 806.0602, 810.5831, 815.106, 819.6289, 824.1518, 828.6747, 833.1976, 837.7205, 842.2434, 846.7663, 851.2892, 855.8121, 860.3564}; 
   
   TH1D *tmpMC_0005_4962 = new TH1D("tmpMC_0005_4962","Other Target",92, xAxis30);
   tmpMC_0005_4962->SetBinContent(7,5.412083);
   tmpMC_0005_4962->SetBinContent(8,2.854513);
   tmpMC_0005_4962->SetBinContent(9,3.988622);
   tmpMC_0005_4962->SetBinContent(10,0.7407119);
   tmpMC_0005_4962->SetBinContent(11,12.62097);
   tmpMC_0005_4962->SetBinContent(12,2.831497);
   tmpMC_0005_4962->SetBinContent(13,4.546255);
   tmpMC_0005_4962->SetBinContent(14,9.24507);
   tmpMC_0005_4962->SetBinContent(15,17.18615);
   tmpMC_0005_4962->SetBinError(7,1.398722);
   tmpMC_0005_4962->SetBinError(8,1.010255);
   tmpMC_0005_4962->SetBinError(9,1.204079);
   tmpMC_0005_4962->SetBinError(10,0.524135);
   tmpMC_0005_4962->SetBinError(11,2.135405);
   tmpMC_0005_4962->SetBinError(12,1.001729);
   tmpMC_0005_4962->SetBinError(13,1.322914);
   tmpMC_0005_4962->SetBinError(14,1.856676);
   tmpMC_0005_4962->SetBinError(15,2.482675);
   tmpMC_0005_4962->SetEntries(164);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetFillColor(ci);
   tmpMC_0005_4962->SetFillStyle(3001);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetLineColor(ci);
   tmpMC_0005_4962->SetLineWidth(5);

   ci = TColor::GetColor("#6600cc");
   tmpMC_0005_4962->SetMarkerColor(ci);
   tmpMC_0005_4962->GetXaxis()->SetTitle("Vertex Z (cm)");
   tmpMC_0005_4962->GetXaxis()->SetRange(8,14);
   tmpMC_0005_4962->GetXaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetXaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetXaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetXaxis()->SetTitleFont(42);
   tmpMC_0005_4962->GetYaxis()->SetTitle("N Events / 1.7 cm");
   tmpMC_0005_4962->GetYaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetYaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetYaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetYaxis()->SetTitleFont(42);
   tmpMC_0005_4962->GetZaxis()->SetLabelFont(42);
   tmpMC_0005_4962->GetZaxis()->SetLabelSize(0.035);
   tmpMC_0005_4962->GetZaxis()->SetTitleSize(0.035);
   tmpMC_0005_4962->GetZaxis()->SetTitleFont(42);
   hs->Add(tmpMC_0005_4962,"");
   hs->Draw("same axis");
   
   TLegend *leg = new TLegend(0.5706,0.554,0.855,0.89,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(10);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("tmpMC_0005_4962","Other Target","f");

   ci = TColor::GetColor("#6600cc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#6600cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC_0004_4962","Downstream Scint.","f");

   ci = TColor::GetColor("#ff9933");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#ff9933");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC_0003_4962","Upstream Scint.","f");

   ci = TColor::GetColor("#cc00cc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#cc00cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC_0002_4962","Iron","f");

   ci = TColor::GetColor("#cc0000");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#cc0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC_0001_4962","Lead","f");

   ci = TColor::GetColor("#006699");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#006699");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   entry=leg->AddEntry("tmpMC_0000_4962","Carbon","f");

   ci = TColor::GetColor("#660066");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3001);

   ci = TColor::GetColor("#660066");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(5);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(62);
   leg->Draw();
   TLatex *   tex = new TLatex(0.18,0.8925,"MINER#nuA Work In Progress");
tex->SetNDC();
   tex->SetTextAlign(13);

   ci = TColor::GetColor("#cc0000");
   tex->SetTextColor(ci);
   tex->SetTextFont(112);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.5,0.945,"DIS Sample - Iron of Target 2");
tex->SetNDC();
   tex->SetTextAlign(21);
   tex->SetTextSize(0.042);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7,0.67,"");
tex->SetNDC();
   tex->SetTextAlign(22);
   tex->SetLineWidth(2);
   tex->Draw();
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->Modified();
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->cd();
   RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL->SetSelected(RegionPosStackPlot_vtxz_t02_z26_minervame6A_v1_CaffeMENNDL);
}
