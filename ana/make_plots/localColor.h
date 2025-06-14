#include "TColor.h"

//======================================================================
enum EColorPalettes
{
  kAlphabetPalette, kKellyPalette, k36Palette, kGlasbeyPalette, kBrewerSet1Palette, kBrewerDark2Palette, kdarker36Palette
};

//======================================================================
const std::vector<int>& getColors(int whichColours)
{
  const int alphabetColors[]={
    TColor::GetColor(  25,  25,  25 ), //ebony
    TColor::GetColor( 240, 163, 255 ), //amethyst
    TColor::GetColor( 255,   0,  16 ), //red
    TColor::GetColor( 128, 128, 128 ), //iron
    TColor::GetColor(  25, 164,   5 ), //orpiment
    TColor::GetColor(   0, 153, 143 ), //turquoise
    TColor::GetColor(   0,  51, 128 ), //navy
    TColor::GetColor(  94, 241, 242 ), //sky
    TColor::GetColor( 157, 204,   0 ), //lime
    TColor::GetColor( 153,  63,   0 ), //caramel
    TColor::GetColor( 224, 255, 102 ), //uranium
    TColor::GetColor(  76,   0,  92 ), //damson
    TColor::GetColor( 255, 168, 187 ), //pink
    TColor::GetColor( 194,   0, 136 ), //mallow
    TColor::GetColor( 255, 204, 153 ), //honeydew
    TColor::GetColor(  43, 206,  72 ), //green
    TColor::GetColor(   0, 117, 220 ), //blue
    TColor::GetColor(   0,  92,  49 ), //forest
    TColor::GetColor( 255, 225,   0 ), //yellow
    TColor::GetColor( 153,   0,   0 ), //wine
    TColor::GetColor( 143, 124,   0 ), //khaki
    TColor::GetColor(  16,  10, 255 ), //violet
    TColor::GetColor( 255, 255, 128 ), //xanthin
    TColor::GetColor( 255,  80,   0 ), //zinnia
    TColor::GetColor( 148, 255, 181 ), //jade
    TColor::GetColor(  66, 102,   0 )  //quagmire
  };
  
  const int kellyColors[]={
    // TColor::GetColor("#f2f3f4"), // white,
    TColor::GetColor("#222222"), // black,
    TColor::GetColor("#f3c300"), // yellow,
    TColor::GetColor("#875692"), // purple,
    TColor::GetColor("#f38400"), // orange,
    TColor::GetColor("#a1caf1"), // lightblue,
    TColor::GetColor("#be0032"), // red,
    TColor::GetColor("#c2b280"), // buff,
    TColor::GetColor("#848482"), // gray,
    TColor::GetColor("#008856"), // green,
    TColor::GetColor("#e68fac"), // purplishpink,
    TColor::GetColor("#0067a5"), // blue,
    TColor::GetColor("#f99379"), // yellowishpink,
    TColor::GetColor("#604e97"), // violet,
    TColor::GetColor("#f6a600"), // orangeyellow,
    TColor::GetColor("#b3446c"), // purplishred,
    TColor::GetColor("#dcd300"), // greenishyellow,
    TColor::GetColor("#882d17"), // reddishbrown,
    TColor::GetColor("#8db600"), // yellowgreen,
    TColor::GetColor("#654522"), // yellowishbrown,
    TColor::GetColor("#e25822"), // reddishorange,
    TColor::GetColor("#2b3d26") // olivegreen
  };

  const int palette36Colors[]={
    TColor::GetColor("#5A5156"),//0 grey
    // TColor::GetColor("#E4E1E3"), // much too light
    TColor::GetColor("#F6222E"),//1 red
    TColor::GetColor("#FE00FA"),//2 hot pink
    TColor::GetColor("#0ABD24"),//3 lime green    //TColor::GetColor("#16FF32"),//3 lime green    
    TColor::GetColor("#3283FE"),//4 medium blue
    TColor::GetColor("#FEAF16"),//5 amber
    TColor::GetColor("#1C8356"),//6 better dark teal green    
    TColor::GetColor("#0BEDBC"),//7 turqoise //"#1CFFCE"),//7 turqoise, kinda light
    TColor::GetColor("#FF28A7"),//#E51990"),//8 dark pink //TColor::GetColor("#B00068"),//light slighly purpley maroon
    TColor::GetColor("#2ED9FF"),//9 light blue
    TColor::GetColor("#AA0DFE"),//10(old 11) purple     
    TColor::GetColor("#DEA0FD"),//11 slighly pink light lavendar
    //switching the pink to 11 and purple to 10 because I like it better
    //TColor::GetColor("#AA0DFE"),//11 purple
    TColor::GetColor("#F8A19F"),//12 peach
    TColor::GetColor("#325A9B"),//13 light navy
    TColor::GetColor("#C4451C"),//14 burt sienna 
    TColor::GetColor("#1CBE4F"),//15 kelly green
    TColor::GetColor("#85660D"),//ugly greenish light brown
    TColor::GetColor("#B10DA1"),//redish purple
    TColor::GetColor("#FBE426"),//light yellow
    TColor::GetColor("#B5EFB5"),//very light green
    TColor::GetColor("#16FF32"),//lime green //TColor::GetColor("#FA0087"),//red pink
    TColor::GetColor("#FC1CBF"),//purple pink
    TColor::GetColor("#F7E1A0"),//beige
    TColor::GetColor("#C075A6"),//violetish lavender
    TColor::GetColor("#782AB6"),//24 really nice purple violet
    TColor::GetColor("#AAF400"),//chartreus
    TColor::GetColor("#BDCDFF"),//very  pale blue
    TColor::GetColor("#822E1C"),//brownish red
    TColor::GetColor("#90AD1C"),//terrible green color, moving toward end of list
    TColor::GetColor("#7ED7D1"),//light baby blue
    TColor::GetColor("#1C7F93"),//dark sky blue
    TColor::GetColor("#D85FF7"),//lilac
    TColor::GetColor("#683B79"),//dark plum (grayish tint)
    TColor::GetColor("#66B0FF"),//sky blue (without any turquise tints)
    TColor::GetColor("#3B00FB"),//blue-indigo blue
  };
  
  const int darker36Colors[]={    TColor::GetColor("#5A5156"),//0 grey
    // TColor::GetColor("#E4E1E3"), // much too light
    TColor::GetColor("#F6222E"),//1 red
    TColor::GetColor("#FE00FA"),//2 hot pink
    TColor::GetColor("#17C62C"),//3 lime green    //TColor::GetColor("#16FF32"),//3 lime green    
    TColor::GetColor("#3283FE"),//4 medium blue
    TColor::GetColor("#FEAF16"),//5 amber
    TColor::GetColor("#1aaf6e"),//6 better dark teal green    
    TColor::GetColor("#31c6c3"),//"#16e2de"),//7 turqoise //"#1CFFCE"),//7 turqoise, too light
				  //TColor::GetColor("#E51990"),//8 dark pink //
    TColor::GetColor("#B00068"),//8 light slighly purpley maroon
    TColor::GetColor("#2ED9FF"),//9 light blue
    TColor::GetColor("#AA0DFE"),//10(old 11) purple     
    TColor::GetColor("#DEA0FD"),//11 slighly pink light lavendar
    //switching the pink to 11 and purple to 10 because I like it better
    //TColor::GetColor("#AA0DFE"),//11 purple
    TColor::GetColor("#F8A19F"),//12 peach
    TColor::GetColor("#325A9B"),//13 light navy
    TColor::GetColor("#C4451C"),//14 burt sienna 
    TColor::GetColor("#1CBE4F"),//15 kelly green
    TColor::GetColor("#85660D"),//ugly greenish light brown
    TColor::GetColor("#B10DA1"),//redish purple
    TColor::GetColor("#FBE426"),//light yellow
    TColor::GetColor("#B5EFB5"),//very light green
    TColor::GetColor("#16FF32"),//lime green //TColor::GetColor("#FA0087"),//red pink
    TColor::GetColor("#FC1CBF"),//purple pink
    TColor::GetColor("#F7E1A0"),//beige
    TColor::GetColor("#C075A6"),//violetish lavender
    TColor::GetColor("#782AB6"),//24 really nice purple violet
    TColor::GetColor("#AAF400"),//chartreus
    TColor::GetColor("#BDCDFF"),//very  pale blue
    TColor::GetColor("#822E1C"),//brownish red
    TColor::GetColor("#90AD1C"),//terrible green color, moving toward end of list
    TColor::GetColor("#7ED7D1"),//light baby blue
    TColor::GetColor("#1C7F93"),//dark sky blue
    TColor::GetColor("#D85FF7"),//lilac
    TColor::GetColor("#683B79"),//dark plum (grayish tint)
    TColor::GetColor("#66B0FF"),//sky blue (without any turquise tints)
    TColor::GetColor("#3B00FB"),//blue-indigo blue
  };


  const int glasbeyColors[]={
    TColor::GetColor("#0000FF"),
    TColor::GetColor("#FF0000"),
    TColor::GetColor("#00FF00"),
    TColor::GetColor("#000033"),
    TColor::GetColor("#FF00B6"),
    TColor::GetColor("#005300"),
    TColor::GetColor("#FFD300"),
    TColor::GetColor("#009FFF"),
    TColor::GetColor("#9A4D42"),
    TColor::GetColor("#00FFBE"),
    TColor::GetColor("#783FC1"),
    TColor::GetColor("#1F9698"),
    TColor::GetColor("#FFACFD"),
    TColor::GetColor("#B1CC71"),
    TColor::GetColor("#F1085C"),
    TColor::GetColor("#FE8F42"),
    TColor::GetColor("#DD00FF"),
    TColor::GetColor("#201A01"),
    TColor::GetColor("#720055"),
    TColor::GetColor("#766C95"),
    TColor::GetColor("#02AD24"),
    TColor::GetColor("#C8FF00"),
    TColor::GetColor("#886C00"),
    TColor::GetColor("#FFB79F"),
    TColor::GetColor("#858567"),
    TColor::GetColor("#A10300"),
    TColor::GetColor("#14F9FF"),
    TColor::GetColor("#00479E"),
    TColor::GetColor("#DC5E93"),
    TColor::GetColor("#93D4FF"),
    TColor::GetColor("#004CFF")
  };

  const int brewerSet1Colors[]={
    TColor::GetColor("#e41a1c"),
    TColor::GetColor("#377eb8"),
    TColor::GetColor("#4daf4a"),
    TColor::GetColor("#984ea3"),
    TColor::GetColor("#ff7f00"),
    TColor::GetColor("#ffff33"),
    TColor::GetColor("#a65628"),
    TColor::GetColor("#f781bf")
  };

  const int brewerDark2Colors[]={
    TColor::GetColor( 27, 158, 119),
    TColor::GetColor(217,  95,   2),
    TColor::GetColor(117, 112, 179),
    TColor::GetColor(231,  41, 138),
    TColor::GetColor(102, 166,  30),
    TColor::GetColor(230, 171,   2),
    TColor::GetColor(166, 118,  29),
    TColor::GetColor(102, 102, 102)
  };

  const int nPalettes=7;//6;
  static std::vector<std::vector<int> > allPalettes(nPalettes);
  static bool firstCall=true;

  if(firstCall){
    for(int j=0; j<nPalettes; ++j){
      int nColors=0;
      const int* firstColor=NULL;

      switch(j){
      case 0:
        nColors=26;
        firstColor=alphabetColors;
        break;
      case 1:
        nColors=21;
        firstColor=kellyColors;
        break;
      case 2:
        nColors=35; // 36 because I commented one out
        firstColor=palette36Colors;
        break;
      case 6:
        nColors=35; // 36 because I commented one out
        firstColor=darker36Colors;
        break;
      case 3:
        nColors=30;
        firstColor=glasbeyColors;
        break;
      case 4:
        nColors=8;
        firstColor=brewerSet1Colors;
        break;
      case 5:
        nColors=8;
        firstColor=brewerDark2Colors;
        break;
      default:
        assert(0 && "whichColours must be 0-5");
      }

      for(int i=0; i<nColors; ++i){
        allPalettes[j].push_back(firstColor[i]);
      }
    }
    firstCall=false;
  }
  return allPalettes[whichColours];
}

//======================================================================
/* void autoColorHists(TPad* pad, int whichColours=kBrewerSet1Palette) */
/* { */
/*   const std::vector<int>& colours=getColors(whichColours); */

/*   std::vector<TH1*> hists=getPadHists(pad); */
/*   for(unsigned int i=0; i<hists.size(); ++i){ */
/*     hists[i]->SetLineColor(colours[i%colours.size()]); */
/*   } */
/*   pad->Draw(); */
/* } */
