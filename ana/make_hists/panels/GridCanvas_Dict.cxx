// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GridCanvas_Dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "GridCanvas.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_GridCanvas(void *p = nullptr);
   static void *newArray_GridCanvas(Long_t size, void *p);
   static void delete_GridCanvas(void *p);
   static void deleteArray_GridCanvas(void *p);
   static void destruct_GridCanvas(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GridCanvas*)
   {
      ::GridCanvas *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GridCanvas >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("GridCanvas", ::GridCanvas::Class_Version(), "GridCanvas.h", 13,
                  typeid(::GridCanvas), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GridCanvas::Dictionary, isa_proxy, 4,
                  sizeof(::GridCanvas) );
      instance.SetNew(&new_GridCanvas);
      instance.SetNewArray(&newArray_GridCanvas);
      instance.SetDelete(&delete_GridCanvas);
      instance.SetDeleteArray(&deleteArray_GridCanvas);
      instance.SetDestructor(&destruct_GridCanvas);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GridCanvas*)
   {
      return GenerateInitInstanceLocal(static_cast<::GridCanvas*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::GridCanvas*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr GridCanvas::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *GridCanvas::Class_Name()
{
   return "GridCanvas";
}

//______________________________________________________________________________
const char *GridCanvas::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int GridCanvas::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GridCanvas::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GridCanvas::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GridCanvas*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void GridCanvas::Streamer(TBuffer &R__b)
{
   // Stream an object of class GridCanvas.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GridCanvas::Class(),this);
   } else {
      R__b.WriteClassBuffer(GridCanvas::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GridCanvas(void *p) {
      return  p ? new(p) ::GridCanvas : new ::GridCanvas;
   }
   static void *newArray_GridCanvas(Long_t nElements, void *p) {
      return p ? new(p) ::GridCanvas[nElements] : new ::GridCanvas[nElements];
   }
   // Wrapper around operator delete
   static void delete_GridCanvas(void *p) {
      delete (static_cast<::GridCanvas*>(p));
   }
   static void deleteArray_GridCanvas(void *p) {
      delete [] (static_cast<::GridCanvas*>(p));
   }
   static void destruct_GridCanvas(void *p) {
      typedef ::GridCanvas current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::GridCanvas

namespace {
  void TriggerDictionaryInitialization_GridCanvas_Dict_Impl() {
    static const char* headers[] = {
"GridCanvas.h",
nullptr
    };
    static const char* includePaths[] = {
"/exp/minerva/app/users/omorenop/cmtuser/git-Mat/MINERvA101/opt/lib/../include/",
"/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/include/root",
"/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/ana/make_hists/panels/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "GridCanvas_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$GridCanvas.h")))  GridCanvas;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "GridCanvas_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "GridCanvas.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"GridCanvas", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("GridCanvas_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_GridCanvas_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_GridCanvas_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_GridCanvas_Dict() {
  TriggerDictionaryInitialization_GridCanvas_Dict_Impl();
}
