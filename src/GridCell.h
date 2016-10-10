//*****************************************************************************
// TGridCell . H :
// written by Tom Harwood Jan-Aug 2011
// base class for data which is stored for land.. (NULL for sea)
// Handles multiple time steps via juggled pairs of stacks, some of which are
// dynamically allocated.
// Use of increment to interpolate between known data steps
// Can take allocated species data
//
//*****************************************************************************
#ifndef GridCellH
#define GridCellH

#include <fstream>

using std::ifstream;

class TGridCell
  {
  public:
   // constructor
   TGridCell(void);
   // destructor
   ~TGridCell(void);
   //GET
   /*
   float* GetGDMData(const int arg_i_stack){return f_gdm_data_ptr[arg_i_stack];}
   float GetCondition(const int arg_i_which){return (f_condition[arg_i_which]);}  // arg_i_which  From, To      ( 0 or 1)
   short int GetRichness(const int arg_i_which){return ((short int)f_richness[arg_i_which]);}  // rounds down all the way..we must have whole species!
   int IsFrozen(void){if (f_gdm_data_ptr[K_FREEZE][0]!=K_MELT)return TRUE; else return FALSE;}
   unsigned short int GetInflux(void){return si_influx;}
   unsigned short int GetOutflux(void){return si_outflux;}
   short int* GetSpeciesData(void){return si_spp_data;}
   short int* GetKnownSpeciesData(void){return si_known_spp_data;} // KM Oct 2012 //
   */
   float* GetEnvData(void){return f_env_data;} // KM May 2014 //
   double* GetCurrentZ(void){return d_env_current_z;} // KM May 2014 //
   double* GetCurrentVp(void){return d_env_current_Vp;} // KM May 2014 //
   double* GetNextZ(void){return d_env_next_z;} // KM May 2014 //
   double* GetNextVp(void){return d_env_next_Vp;} // KM May 2014 //
   short int GetOccurrence(void){return s_present;}   // KM May 2014 //
   short int GetNewlyColonised(void){return s_newly_colonised;}   // KM May 2014 //
   short int GetTimeLastPresent(void){return s_time_last_present;} // KM May 2014 //

   //SET
   void SetEnvData(float* arg_f_data_ptr){f_env_data=arg_f_data_ptr;}    // KM May 2014 //
   void SetCurrentZ(double* arg_d_data_ptr){d_env_current_z=arg_d_data_ptr;} // KM May 2014 //
   void SetCurrentVp(double* arg_d_data_ptr){d_env_current_Vp=arg_d_data_ptr;}  // KM May 2014 //
   void SetNextZ(double* arg_d_data_ptr){d_env_next_z=arg_d_data_ptr;} // KM May 2014 //
   void SetNextVp(double* arg_d_data_ptr){d_env_next_Vp=arg_d_data_ptr;}  // KM May 2014 //
   void SetOccurrence(short int arg_s_value){s_present=arg_s_value;}            // KM June 2014 //
   void SetNewlyColonised(short int arg_s_value){s_newly_colonised=arg_s_value;}            // KM June 2014 //
   void SetTimeLastPresent(short int arg_s_value){s_time_last_present=arg_s_value;}            // KM June 2014 //
   void SetENVStack(ifstream* arg_InFileStream_ptr,
					const int arg_n_layers);

   // For Adaptor...
   float* f_env_data;  		// pointer to a dynamically allocated array to hold the values of each environmental variable for this cell
   double* d_env_current_z; // pointer to a dynamically allocated array to hold the values of the mean tolerance attribute (z) for each environment variable
   double* d_env_current_Vp; // pointer to a dynamically allocated array to hold the values of the variance in tolerance attribute (Vp) for each environment variable
   double* d_env_next_z; // pointer to a dynamically allocated array to hold the values of the mean tolerance attribute (z) for each environment variable at the next time point
   double* d_env_next_Vp; // pointer to a dynamically allocated array to hold the values of the variance in tolerance attribute (Vp) for each environment variable at the next time point
   short int s_time_last_present; // holds the last timestep this species was present in this grid cell
   short int s_present;         // holds whether the species is present (1) or absent (0) in this gris cell
   short int s_newly_colonised;

  protected:
  }; // end class TGridCell
typedef TGridCell* TGridCell_ptr;
//---------------------------------------------------------------------------
#endif
