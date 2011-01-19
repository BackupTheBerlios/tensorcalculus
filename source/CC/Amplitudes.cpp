/*
 * Copyright (C) 2011 Henry Auer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef AMPLITUDES_CPP_
#define AMPLITUDES_CPP_

#include "CC/Amplitudes.hpp"
#include "CC/C4Interface.hpp"
#include "Vector/VectorOperators.hpp"

namespace TensorCalculus{

   CPTensorRepresentation<double> readEpsilonInvDKTS (){
      CPTensorRepresentation<double> EpsilonInv("compute/epsilonInv.dkts");

      int occ, virt;
      readOccN(occ, occ, virt, virt);

      int start, end, length;
      if(virt<occ){
         start = 0;
         end = 2;
         length = virt;
      }else if (occ<virt){
         start = 2;
         end = 4;
         length = occ;
      }else{
         return (EpsilonInv);
      }
      std::vector< std::vector <double> > v(EpsilonInv.getV());
      for(int mu = start; mu<end; mu++){
         v[mu].resize(length * EpsilonInv.getSummation(0));
         for(int j=0; j<EpsilonInv.getSummation(0); j++){
            for(int i=0; i<length; i++){
               v[mu][i+length*j] = EpsilonInv.getV()[mu][i+std::max(occ,virt)*j];
            }
         }
      }
      std::vector<int> compDim(4);
      compDim[0] = virt;
      compDim[1] = virt;
      compDim[2] = occ;
      compDim[3] = occ;
      return(CPTensorRepresentation<double> (EpsilonInv.getSummation(0), compDim, v) );
   }

   void V2T (MPSRepresentation<double> &V_mnop){
      CPTensorRepresentation<double> EpsilonInv (readEpsilonInvDKTS());
      V_mnop.setHadamardProduct(EpsilonInv);
      //int occ, virt;
      //readOccN (occ, occ, virt, virt);
      //todo ... occ virt mitDimensions vergleichen
      // --> chem oder phys oder andere notation !
   }

   void Amplitudes::init (const std::vector<double> &spinAdaptV, const std::vector<double> &v){
      readOccN(occ, occ, virt, virt);
      readOrbitalEnergies(f);
      aoCount = occ + virt;
      int size = aoCount * aoCount * aoCount * aoCount;
      fTensor.resize(size);
      size = occ * occ * virt * virt;
      f_vv_oo.resize(size);
      t_vv_oo.resize(size);
      v_vv_oo = spinAdaptV;
      v_vv_oo_simple = v;
   }

   void Amplitudes::makeOrbitalEnergieTensor (){

      int index;
      for(int i=0; i<aoCount; i++){
         for(int j=0; j<aoCount; j++){
            for(int k=0; k<aoCount; k++){
               for(int l=0; l<aoCount; l++){
                  index = l + aoCount * (k + aoCount * (j + aoCount * i));
                  fTensor[index] = -1.0 * f[l] - f[k] + f[i] +f[j];
                  fTensor[index] = 1.0 / fTensor[index];
               }
            }
         }
      }
   }

   void Amplitudes::makeF_vv_oo (){
      int index;
     for(int i=0; i<occ; i++){
        for(int j=0; j<occ; j++){
           for(int k=0; k<virt; k++){
              for(int l=0; l<virt; l++){
                 index = l + virt * (k + virt * (j + occ * i));
                 f_vv_oo[index] = -1.0 * f[l+occ] - f[k+occ] + f[j] +f[i];
                 f_vv_oo[index] = 1.0 / f_vv_oo[index];
              }
           }
        }
     }
   }

   void Amplitudes::makeT_vv_oo (){
      int size = occ * occ * virt * virt;
      for(int i=0; i<size; i++){
         t_vv_oo[i] = v_vv_oo[i] * f_vv_oo[i];
         //std::cout<<t_vv_oo[i]<<"  ";
      }
      //std::cout<<std::endl;
   } //end makeT ()

   double Amplitudes::MP2Energie (){
      //(*this).makeOrbitalEnergieTensor();
      (*this).makeF_vv_oo();
      (*this).makeT_vv_oo();
      int size = occ*occ*virt*virt;
      double temp = 0.0;
      for(int i=0; i<size; i++){
         temp += v_vv_oo_simple[i]*t_vv_oo[i];
      }
      //temp /= 2.0;
      return(temp);
      //return(Blas<double>::dot(size, &v_vv_oo[0], 1, &t_vv_oo[0], 1));
   }//end MP2Energie
/*
   void Amplitudes::makeOrbitalEnergieCPTensor (){
      // code von Udo, modifiziert
      int virt, occ;
      readOccN(occ, occ, virt, virt);
      int aoCount = virt + occ;

      std::vector<double> componentDimensions(4);
      componentDimensions[0] = aoCount;
      componentDimensions[1] = aoCount;
      componentDimensions[2] = aoCount;
      componentDimensions[3] = aoCount;

      int r = 42;

      std::vector<double> alpha(r);
      std::vector<double> omega(r);

      alpha[0]  =  0.0000000001504533796769452920001572409241;
      alpha[1]  =  0.0000000011531074953405321307662126075555;
      alpha[2]  =  0.0000000048888187124083736735983509053166;
      alpha[3]  =  0.0000000164913263438988012605406899407238;
      alpha[4]  =  0.0000000482281745320224058462402669281622;
      alpha[5]  =  0.0000001273504255333853829028322031240319;
      alpha[6]  =  0.0000003110894054663478634846938080616295;
      alpha[7]  =  0.0000007143105683543448472755047572440857;
      alpha[8]  =  0.0000015590485153608081304051440251795698;
      alpha[9]  =  0.0000032610487264231998748898275792626879;
      alpha[10] =  0.0000065777069965383905195249778677341968;
      alpha[11] =  0.0000128561535096668727541317751632222923;
      alpha[12] =  0.0000244422071323536523197078529224114500;
      alpha[13] =  0.0000453442101154741060146244694268853961;
      alpha[14] =  0.0000822956578083203576276421554493156751;
      alpha[15] =  0.0001464346732174451759210819678418483147;
      alpha[16] =  0.0002559275716031920948236231508517157796;
      alpha[17] =  0.0004400233198752026714698137247871317923;
      alpha[18] =  0.0007452571228546059473083139855758422199;
      alpha[19] =  0.0012448547035166384835101839102279586147;
      alpha[20] =  0.0020528658296695951231062817398873021624;
      alpha[21] =  0.0033452339853937848513022501649796791590;
      alpha[22] =  0.0053909681941535552222664179605055023714;
      alpha[23] =  0.0085979316261365772021513470373121901247;
      alpha[24] =  0.0135796481070218221221172316508229993559;
      alpha[25] =  0.0212521535780708742816496941915871410345;
      alpha[26] =  0.0329735581107251665781412638955849558897;
      alpha[27] =  0.0507440050739478467879900805415083198113;
      alpha[28] =  0.0774906195750688788384898976480119614507;
      alpha[29] =  0.1174715178598124021522330386047006101080;
      alpha[30] =  0.1768459776869590412189331765646649330392;
      alpha[31] =  0.2644759110627820372706013596353358252600;
      alpha[32] =  0.3930492819191733021961469157767865567621;
      alpha[33] =  0.5806538052951022052224419600641169836308;
      alpha[34] =  0.8529897882443225512580538272278118938630;
      alpha[35] =  1.2465211951108858604280599235991644491151;
      alpha[36] =  1.8130974497943952231104880468137707794085;
      alpha[37] =  2.6271499107685698829687398481169680053426;
      alpha[38] =  3.7981600722550915947162836205919234089379;
      alpha[39] =  5.4962274554572003331699703299051407157094;
      alpha[40] =  8.0191886922357356179066889545481444656616;
      alpha[41] = 12.0573866304389756450021753408918812056072;

      omega[0]  =   0.0000000004219320915572596232365877921406;
      omega[1]  =   0.0000000018684214423633261599764249918123;
      omega[2]  =   0.0000000063980839195715625748534346181758;
      omega[3]  =   0.0000000187278992702136833688625133023618;
      omega[4]  =   0.0000000490502922225399008082726599651067;
      omega[5]  =   0.0000001182929838000904322672873366452947;
      omega[6]  =   0.0000002675383524852711308788151492505010;
      omega[7]  =   0.0000005745234884851597950263092403360061;
      omega[8]  =   0.0000011818669725721989293418032504373355;
      omega[9]  =   0.0000023443715607824973478721145007765202;
      omega[10] =   0.0000045068681726955300712347635451595707;
      omega[11] =   0.0000084302640191914659673003302782902466;
      omega[12] =   0.0000153926953424398783918871521538668766;
      omega[13] =   0.0000275063282148737354199025115296210753;
      omega[14] =   0.0000482103260726183562083977135091176353;
      omega[15] =   0.0000830295200990567342604392725414173826;
      omega[16] =   0.0001407301592452955568341578016888608615;
      omega[17] =   0.0002350640278093540352986928182303129342;
      omega[18] =   0.0003873774517226005101409423302030038405;
      omega[19] =   0.0006304822293962166182500015266642473533;
      omega[20] =   0.0010143549460067730317708843430085652315;
      omega[21] =   0.0016144680116281527549437402779992922763;
      omega[22] =   0.0025438852570893892895838138282563201109;
      omega[23] =   0.0039707110037390608483969596693985426583;
      omega[24] =   0.0061431099574655103204423177605530970169;
      omega[25] =   0.0094249775610127827048373145950310725283;
      omega[26] =   0.0143465193386527931298203499162202678718;
      omega[27] =   0.0216756054977696673311476178655499147396;
      omega[28] =   0.0325179589832321061038458709369081134355;
      omega[29] =   0.0484572373295624925637473719847170272601;
      omega[30] =   0.0717502414883026883528870099038243779432;
      omega[31] =   0.1055984902416372645539703774286710569186;
      omega[32] =   0.1545266611196524010634527693253126301443;
      omega[33] =   0.2249144111603383853754782595313344017995;
      omega[34] =   0.3257601636970398329936401266015977284951;
      omega[35] =   0.4698296803142691221901600107102581205254;
      omega[36] =   0.6755350077022725790174342841432064687979;
      omega[37] =   0.9704388241513622397017499454641153988632;
      omega[38] =   1.3990172699260266954814427298181556125201;
      omega[39] =   2.0437808601951682316213293466589107083564;
      omega[40] =   3.1008805126561835041089365549993317472399;
      omega[41] =   5.3312998445317940289685243460127139769611;

      //todo ....

   } //end makeOrbitalEnergieCPTensor ()
   */
   CPTensorRepresentation<double> makeEpsilonInv(const std::vector<char> &occState){
      // nach Udo Benedikt
      int occCount, virtCount;
      readOccN(occCount, occCount, virtCount, virtCount);

      int r = 42;

      std::vector< std:: vector <double> > v(4);
      std::vector<int> componentDimensions(4);
      for(int mu=0; mu<4; mu++){
         if (occState[mu] == 'v'){
            v[mu].resize(virtCount * r);
            componentDimensions[mu] = virtCount;
         } else if (occState[mu] == 'o'){
            v[mu].resize(occCount * r);
            componentDimensions[mu] = occCount;
         } else {
            throw std::invalid_argument("undefined occumpation State");
         }
      }

      std::vector<double> alpha(r);
      std::vector<double> omega(r);

      alpha[0]  =  0.0000000001504533796769452920001572409241;
      alpha[1]  =  0.0000000011531074953405321307662126075555;
      alpha[2]  =  0.0000000048888187124083736735983509053166;
      alpha[3]  =  0.0000000164913263438988012605406899407238;
      alpha[4]  =  0.0000000482281745320224058462402669281622;
      alpha[5]  =  0.0000001273504255333853829028322031240319;
      alpha[6]  =  0.0000003110894054663478634846938080616295;
      alpha[7]  =  0.0000007143105683543448472755047572440857;
      alpha[8]  =  0.0000015590485153608081304051440251795698;
      alpha[9]  =  0.0000032610487264231998748898275792626879;
      alpha[10] =  0.0000065777069965383905195249778677341968;
      alpha[11] =  0.0000128561535096668727541317751632222923;
      alpha[12] =  0.0000244422071323536523197078529224114500;
      alpha[13] =  0.0000453442101154741060146244694268853961;
      alpha[14] =  0.0000822956578083203576276421554493156751;
      alpha[15] =  0.0001464346732174451759210819678418483147;
      alpha[16] =  0.0002559275716031920948236231508517157796;
      alpha[17] =  0.0004400233198752026714698137247871317923;
      alpha[18] =  0.0007452571228546059473083139855758422199;
      alpha[19] =  0.0012448547035166384835101839102279586147;
      alpha[20] =  0.0020528658296695951231062817398873021624;
      alpha[21] =  0.0033452339853937848513022501649796791590;
      alpha[22] =  0.0053909681941535552222664179605055023714;
      alpha[23] =  0.0085979316261365772021513470373121901247;
      alpha[24] =  0.0135796481070218221221172316508229993559;
      alpha[25] =  0.0212521535780708742816496941915871410345;
      alpha[26] =  0.0329735581107251665781412638955849558897;
      alpha[27] =  0.0507440050739478467879900805415083198113;
      alpha[28] =  0.0774906195750688788384898976480119614507;
      alpha[29] =  0.1174715178598124021522330386047006101080;
      alpha[30] =  0.1768459776869590412189331765646649330392;
      alpha[31] =  0.2644759110627820372706013596353358252600;
      alpha[32] =  0.3930492819191733021961469157767865567621;
      alpha[33] =  0.5806538052951022052224419600641169836308;
      alpha[34] =  0.8529897882443225512580538272278118938630;
      alpha[35] =  1.2465211951108858604280599235991644491151;
      alpha[36] =  1.8130974497943952231104880468137707794085;
      alpha[37] =  2.6271499107685698829687398481169680053426;
      alpha[38] =  3.7981600722550915947162836205919234089379;
      alpha[39] =  5.4962274554572003331699703299051407157094;
      alpha[40] =  8.0191886922357356179066889545481444656616;
      alpha[41] = 12.0573866304389756450021753408918812056072;

      omega[0]  =   0.0000000004219320915572596232365877921406;
      omega[1]  =   0.0000000018684214423633261599764249918123;
      omega[2]  =   0.0000000063980839195715625748534346181758;
      omega[3]  =   0.0000000187278992702136833688625133023618;
      omega[4]  =   0.0000000490502922225399008082726599651067;
      omega[5]  =   0.0000001182929838000904322672873366452947;
      omega[6]  =   0.0000002675383524852711308788151492505010;
      omega[7]  =   0.0000005745234884851597950263092403360061;
      omega[8]  =   0.0000011818669725721989293418032504373355;
      omega[9]  =   0.0000023443715607824973478721145007765202;
      omega[10] =   0.0000045068681726955300712347635451595707;
      omega[11] =   0.0000084302640191914659673003302782902466;
      omega[12] =   0.0000153926953424398783918871521538668766;
      omega[13] =   0.0000275063282148737354199025115296210753;
      omega[14] =   0.0000482103260726183562083977135091176353;
      omega[15] =   0.0000830295200990567342604392725414173826;
      omega[16] =   0.0001407301592452955568341578016888608615;
      omega[17] =   0.0002350640278093540352986928182303129342;
      omega[18] =   0.0003873774517226005101409423302030038405;
      omega[19] =   0.0006304822293962166182500015266642473533;
      omega[20] =   0.0010143549460067730317708843430085652315;
      omega[21] =   0.0016144680116281527549437402779992922763;
      omega[23] =   0.0039707110037390608483969596693985426583;
      omega[22] =   0.0025438852570893892895838138282563201109;
      omega[24] =   0.0061431099574655103204423177605530970169;
      omega[25] =   0.0094249775610127827048373145950310725283;
      omega[26] =   0.0143465193386527931298203499162202678718;
      omega[27] =   0.0216756054977696673311476178655499147396;
      omega[28] =   0.0325179589832321061038458709369081134355;
      omega[29] =   0.0484572373295624925637473719847170272601;
      omega[30] =   0.0717502414883026883528870099038243779432;
      omega[31] =   0.1055984902416372645539703774286710569186;
      omega[32] =   0.1545266611196524010634527693253126301443;
      omega[33] =   0.2249144111603383853754782595313344017995;
      omega[34] =   0.3257601636970398329936401266015977284951;
      omega[35] =   0.4698296803142691221901600107102581205254;
      omega[36] =   0.6755350077022725790174342841432064687979;
      omega[37] =   0.9704388241513622397017499454641153988632;
      omega[38] =   1.3990172699260266954814427298181556125201;
      omega[39] =   2.0437808601951682316213293466589107083564;
      omega[40] =   3.1008805126561835041089365549993317472399;
      omega[41] =   5.3312998445317940289685243460127139769611;

      std::vector<double> orbitalEnergies (occCount + virtCount);
      readOrbitalEnergies (orbitalEnergies);

      using namespace VectorOperators;

      double HOMO = orbitalEnergies[occCount - 1];
      double LUMO = orbitalEnergies[occCount];
      double epsilon_min = 2.0*(HOMO + LUMO);

      for(int mu = 0; mu < 4; mu++){
         for (int j=0; j < r; j++){
            if(occState[mu] == 'o'){
               for(int i=0; i < occCount; i++){
                  if(mu == 0){
                     v[mu][i + occCount * j] = -1.0 * omega[j] / epsilon_min
                           * exp(-alpha[j] / epsilon_min * orbitalEnergies[i]);
                  } else {
                     v[mu][i + occCount * j] = exp(-alpha[j] / epsilon_min * orbitalEnergies[i]);
                  }
               }
            } else {
               for(int i=0; i < virtCount; i++){
                  if(mu == 0){
                     v[mu][i + virtCount * j] = -1.0 * omega[j] / epsilon_min
                           * exp(-alpha[j] / epsilon_min * orbitalEnergies[i+occCount]);
                  } else {
                     v[mu][i + virtCount * j] = exp(-alpha[j] / epsilon_min * orbitalEnergies[i+occCount]);
                  }
               }
            }
         }
      }//end for(mu=0; mu < d; mu++)

      return(CPTensorRepresentation<double> (r, componentDimensions, v));
   }

   void V2T (MPSRepresentation<double> &V_mnop,const std::vector<char> &occState){
      V_mnop.setHadamardProduct(makeEpsilonInv(occState));
   }
}// end namespace TensorCalculus


#endif /* AMPLITUDES_CPP_ */
