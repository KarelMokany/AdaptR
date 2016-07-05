//---------------------------------------------------------------------------

//*****************************************************************************
//*
//* ADAPTOR  FILE UTILS
//*
//* Miscellaneous funcs to do stuff for Adaptor that don't really fit anywhere else
//*
//*****************************************************************************
#ifndef AdaptorFileUtils_V2H
#define AdaptorFileUtils_V2H
using std::string;

class TRaft;

//****************************************************************************
//* LOAD PARAMS FILE

//*  Returns Success or Failure 1 or 0
//***************************************************************************
int Run_Adaptor(string arg_filename);

double phi(double arg_d_x);

double calculate_survival_percentage(double arg_d_x,
									 double arg_d_Z,
									 double arg_d_Vp);

double calculate_survival_percentage_LowLimit(double arg_d_x,
											  double arg_d_Z,
											  double arg_d_Vp);

double calculate_above_limit_proportion(double arg_d_x,
										double arg_d_Z,
										double arg_d_Vp);

double calculate_below_limit_proportion(double arg_d_x,     //x
										double arg_d_Z,     //mean
										double arg_d_Vp);    //variance

void Calculate_upper_threshold_correction(int arg_i_n_segments,
										  int arg_i_index_for_mean,
										  double arg_d_variance_limit,
										  double* arg_d_revised_Z_upper_limit,		// gets filled here
										  double* arg_d_revised_Vp_upper_limit);   // gets filled here

int round_double(double arg_d_number);

void Implement_dispersal_mode_1_correction(int arg_n_rows,
										   int arg_n_cols,
										   int arg_i_time,
										   int arg_i_n_dispersal_categories,
										   int** arg_i_dispersal_category,
										   short** arg_s_time_last_present,
										   double** arg_d_current_z,     // Gets changed here
										   double** arg_d_current_Vp);    // Gets changed here


void Implement_dispersal_mode_2_correction(int arg_n_rows,
										   int arg_n_cols,
										   int arg_n_rel,
										   int arg_i_time,
										   short** arg_s_time_last_present,
										   int** arg_i_rel_row_col,
										   double* arg_d_rel_prob,
										   double** arg_d_current_z,     // Gets changed here
										   double** arg_d_current_Vp);    // Gets changed here

void Assess_Persistence_and_Adaptation(TRaft* Raft,
									   int arg_n_rows,
									   int arg_n_cols,
									   int arg_n_environments,
									   int arg_i_time,
									   int arg_i_index_for_mean,
									   int arg_i_n_segments,
									   double arg_d_spp_survival_min,
									   double arg_d_increment_per_index,
									   int* arg_b_env_low_adaptation,
									   int* arg_b_env_high_adaptation,
									   float* arg_f_env_low_fund_limit,
									   float* arg_f_env_high_fund_limit,
									   double* arg_d_adapt_limit,
									   double* arg_d_revised_Z_upper_limit,
									   double* arg_d_revised_Vp_upper_limit,
									   double* arg_d_spp_H2,
									   double* arg_d_nat_sel_slope);

void Implement_Dispersal_To_Unoccupied_Cells(TRaft* Raft,
											 int arg_n_rows,
											 int arg_n_cols,
											 int arg_n_rel,
											 int arg_n_environments,
											 int* arg_b_env_low_adaptation,
											 int* arg_b_env_high_adaptation,
											 int** arg_i_rel_row_col,
											 float* arg_f_env_low_fund_limit,
											 float* arg_f_env_high_fund_limit,
											 double* arg_d_rel_prob,
											 double* arg_d_adapt_limit);

void Calculate_Newly_Colonised_Z_Vp(TRaft* Raft,
									int arg_i_row,
									int arg_i_col,
									int arg_n_rel,
									int arg_n_environments,
									int arg_i_ncats_env,
									double arg_d_prob_sum,
									int* arg_b_env_low_adaptation,
									int* arg_b_env_high_adaptation,
									double* arg_d_summed_prob,
									double* arg_d_rel_prob,
									double* arg_d_adapt_limit,
									double** arg_d_rel_Z,
									double** arg_d_rel_Vp);

void Implement_Adaptive_Parameter_Averaging(TRaft* Raft,
											int arg_n_rows,
											int arg_n_cols,
											int arg_n_rel,
											int arg_n_environments,
											double arg_d_resident_weighting,
											int* arg_b_env_low_adaptation,
											int* arg_b_env_high_adaptation,
											int** arg_i_rel_row_col,
											double* arg_d_rel_prob,
											double* arg_d_adapt_limit);

void convert_colonised_cellls_to_occupied(TRaft* Raft,
										  int arg_n_rows,
										  int arg_n_cols);

float uniformRandomDeviate(float arg_f_min,
							float arg_f_max);

/*
int GenerateCircularFill(const int arg_i_max_radius,
						  const int arg_i_min_radius,
						  short ** arg_si_rel_grid);
*/

void Calculate_Newly_Colonised_Z_Vp_Mixture(TRaft* Raft,
											int arg_i_row,
											int arg_i_col,
											int arg_n_rel,
											int arg_n_environments,
											int* arg_b_env_low_adaptation,
											int* arg_b_env_high_adaptation,
											float* arg_f_env_low_fund_limit,
											float* arg_f_env_high_fund_limit,
											double* arg_d_rel_prob,
											double* arg_d_adapt_limit,
											double** arg_d_rel_Z,
											double** arg_d_rel_Vp);

void Implement_Adaptive_Parameter_Averaging_Mixture(TRaft* Raft,
													int arg_n_rows,
													int arg_n_cols,
													int arg_n_rel,
													int arg_n_environments,
													double arg_d_resident_weighting,
													int* arg_b_env_low_adaptation,
													int* arg_b_env_high_adaptation,
													int** arg_i_rel_row_col,
													double* arg_d_rel_prob,
													double* arg_d_adapt_limit);

void  Generate_occurrence_summary(TRaft* Raft,
								  int arg_n_rows,
								  int arg_n_cols,
								  int arg_n_environments,
								  int* arg_b_env_low_adaptation,
								  int* arg_b_env_high_adaptation,
								  short** arg_s_initial_occurrence,
								  double* arg_d_mean_Z_env,
								  double* arg_d_mean_Vp_env,
								  int* arg_i_n_occupied_ptr,
								  int* arg_i_n_initial_occupied_ptr);

//---------------------------------------------------------------------------
#endif
