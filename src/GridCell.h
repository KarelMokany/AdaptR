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
   /*
   void SetCondition(const int arg_i_which,     //Assume From or To    ( either 0 or 1, but can swap)
					 const float arg_f_value){f_condition[arg_i_which]=arg_f_value;}
   void SetRichness(const int arg_i_which,     //Assume From or To
					const float arg_f_value){f_richness[arg_i_which]=arg_f_value;}
   void SetSpeciesData(short int* arg_si_data_ptr){si_spp_data=arg_si_data_ptr;}
   void SetKnownSpeciesData(short int* arg_si_data_ptr){si_known_spp_data=arg_si_data_ptr;} // KM Oct 2012 //
   void SetInflux(const unsigned short int arg_i_value){si_influx=arg_i_value;}
   void SetOutflux(const unsigned short int arg_i_value){si_outflux=arg_i_value;}
   void AddToInflux(const unsigned short int arg_i_value){si_influx+=arg_i_value;}
   void AddToOutflux(const unsigned short int arg_i_value){si_outflux+=arg_i_value;}

   void SetFreeze(const float arg_f_freeze,    // K_FREEZE or K_MELT
                  const int arg_i_from_stack,  // can be anything if MELTING otherwise vital
                  const int arg_n_layers);
   void SetGDMStack(const int arg_i_GDMStack,
                    ifstream* arg_InFileStream_ptr,
                    const int arg_n_layers);
   //Allocation
   void LoadGDMStack(const int arg_i_GDMStack,   // Allocates and sets
                     ifstream* arg_InFileStream_ptr,
                     const int arg_n_layers);
   void AllocateGDMStack(const int arg_i_GDMStack,
                         const int arg_n_layers);  // allocates space for the stack at this cell
   void FreeGDMStack(const int arg_i_GDMStack);  // frees space for the stack at this cell
   void AllocateGDMIncrement(const int arg_n_layers);  // allocates space for the stack at this cell
   void FreeGDMIncrement(void); // Free up array of floats
   

   //Incrementationisation
   int CalcGDMIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                        const int arg_n_layers,
                        const int arg_i_from_stack);
   void IncrementGDMIncrement(const int arg_i_to_stack,
                              const int arg_n_layers);
   int CalcCondIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                         const int arg_i_from_stack);
   void IncrementCondIncrement(const int arg_i_to_stack);
   int CalcRichIncrement(const int arg_i_step,   // Calculates arg_i_step increment between each stack
                         const int arg_i_from_stack);
   void IncrementRichIncrement(const int arg_i_to_stack);


   // Calc Similarity overloaded function so we can opt to get ecol dist
   float CalcSimilarity(float* arg_f_OtherData,  // stack to compare with
                        const int arg_i_stack,   //which stack we are
                        const int arg_n_layers,       // how many
                        const float arg_f_intercept);

   float CalcSimilarity(float* arg_f_OtherData,  // stack to compare with
                        const int arg_i_stack,   //which stack we are
                        const int arg_n_layers,  // how many
                        float* arg_f_ecol_dist,   // pointer to float which we can set to the ecol distance
                        const float arg_f_intercept);
      // Calc Dissimilarity overloaded function so we can opt to get ecol dist
   float CalcDissimilarity(float* arg_f_OtherData,  // stack to compare with
                            const int arg_i_stack,   //which stack we are
                            const int arg_n_layers,  // how many
                            const float arg_f_intercept);   // gdm intercept

   float CalcDissimilarity( float* arg_f_OtherData,  // stack to compare with
                            const int arg_i_stack,   //which stack we are
							const int arg_n_layers,  // how many
                            float* arg_f_ecol_dist,   // pointer to float which we can set to the ecol distance
                            const float arg_f_intercept); // gdm intercept

   // keep data public
   float* f_gdm_data_ptr[3]; // array of pointers for the start of dynamically allocated data arrays.
                             // si_gdm_data_ptr[0] From Data, si_gdm_data_ptr[1] ToData, si_data_ptr Condition grids
                             // no point in dynamically allocating..
                             // if a GridCell is non- null we will need all three at some point
                             //  2: frozen from freeze layer to store any backed up data
   float* f_gdm_increment;   // pointer to a dynamically allocated arrayto hold the annual increment (may not be present)
   float f_richness[2];  // data points for richness grids// data points for condition grids
   float f_rich_increment;  // sub step increment to interpolate richness
   float f_condition[2]; // data points for condition grids
   float f_cond_increment;   // sub step increment to interpolate for condition
   short int* si_spp_data;   // pointer to an array of data stored for this occupied cell.
   short int* si_known_spp_data;   // KM Oct 2012 // pointer to an array of known species data stored for this occupied cell. Uses f_richness[1] to indicate length
   unsigned short int si_influx;
   unsigned short int si_outflux;
   */
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
