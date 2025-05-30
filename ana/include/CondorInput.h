#ifndef MNV_CONDORUTILS_CONDORINPUT
#define MNV_CONDORUTILS_CONDORINPUT 1

//forwards
class TString;
class TSystem;

namespace CondorUtils
{
  class CondorInput
  {
    public:
      //! Default constructor
      CondorInput();

      //! Default destructor
      virtual ~CondorInput();

      //! singleton gettor
      static CondorInput& Get();
   
      //! Simple check for the existence of a file
      bool FileExists( const TString& fname ) const;

      //! Are we on the grid?
      bool IsOnGrid() const;

      /*! Where is the directory that holds condor input files?
        @return Expansion of $CONDOR_DIR_INPUT ( returnVal.IsNull() is true if CONDOR_DIR_INPUT is not defined)
       */
      TString GetCondorDirInput() const;

      /*! Where is this condor directory?
        @param[in] 
        @return Expansion of $CONDOR_DIR_<dirname> ( returnVal.IsNull() is true if directory is not defined)
       */
      TString GetCondorDir( const TString& dirname ) const;

      /*! Copy a file from bluearc to the condor input directory
        @param[in] fname Full path and filename of file on bluearc
        @param[in] allowOverwrite If set to true, then copy the file to the condor input area, even if a file with this name is already there.  By default is does nothing if the file is already there.
        @return The return value of the command used to copy the file (0 is success)
       */
      int CopyToInput( const TString& fname, bool allowOverwrite = false ) const;

      /*! What is the full path to this input file?
        @param[in] name Full path and filename of file on bluearc
        @return What filename should you use to open this file?
        */
      TString GetInputFilename( const TString& fname ) const;

      /*! Copy the file if necessary.  Return the filename the user should use.
        @param[in] name Full path and filename of the file on bluarc
        @param[in] allowOverwrite If set to true, then copy the file to the condor input area, even if a file with this name is already there.  By default is does nothing if the file is already there.
        @return What filename should you use to open this file?
        */
      TString FetchInput( const TString& fname, bool allowOverwrite = false  ) const;

    private:
      bool onGrid_;              ///< Are we on a grid machine?
      TString *condorDirInput_;  ///< The value of $CONDOR_DIR_INPUT?
      TSystem *tSystem_;         ///< An "instance" of TSystem


  }; // end class CondorInput

  inline bool CondorInput::IsOnGrid()             const { return onGrid_; }

}//end namespace CondorUtils

#endif //MNV_CONDORUTILS_CONDORINPUT
