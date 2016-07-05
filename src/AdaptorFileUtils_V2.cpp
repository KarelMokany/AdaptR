//**********************************************************
// *
// * ADAPTOR  FILE UTILS
// *
// * Miscellaneous funcs to do stuff for ADAPTOR that don't really fit anywhere else
// *
//*************************************************************
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "DynamicArray.h"
#include "defines.h"
#include "Raft.h"
#include "GridCell.h"
#include "AdaptorFileUtils_V2.h"

using namespace std;
//****************************************************************************
//------------------------------------------------------------------------------
// Run Adaptor is the primary function that reads in the parameter file, performs
// the adaptation calculations, then writes the output files.
int Run_Adaptor(string arg_filename)
{
  // Declare local variables
  int i_row, i_col, i_time, i_loc, i_rel, i_env;
  int n_cols, n_rows;
  int n_timepoints;
  int n_environments;
  int n_spp=1;
  int n_locations;
  int n_rel;
  int i_check = 0;
  int i_index_for_mean;
  int i_buf;
  int i_number;
  int i_data_read;
  int n_occupied;
  int n_initial_occupied;
  int n_not_initial_occupied;
  int i_now;
  int i_runtime;
  int b_write_each_step = 0;
  short int s_present_time;
  short int s_present;
  short int s_new;
  short int s_occurrence;
  short int s_one = 1;
  short int s_zero = 0;
  short int s_nodata = -9999;
  long int i_this_spp;
  long int i_this_env;
  float f_data;
  float f_Env1;
  float f_total_time;
  float f_one = 1;
  float f_nodata = -9999;
  float f_Z;
  double d_data;
  double d_spp_survival_min;
  double d_S;
  double d_i;
  double d_pt_Latitude;
  double d_pt_Longitude;
  double d_yulcorner;
  double d_increment_per_index;
  double d_resident_weighting;
  double d_nodata = -9999;

  // Variables to set for the upper limit correction
  int i_n_segments = 81;
  double d_variance_limit = 4;

  // Allocate ascii file header stuff
  int n_totalcols, n_totalrows;
  double d_xllcorner, d_yllcorner;
  double d_cellsize;
  int i_nodata;

  // Allocate strings
  string str_buffer;
  string str_input_folder_path;
  string str_env1_grids_folder_path;
  string str_spp_grids_folder_path;
  string str_output_folder_path;
  string str_output_rootname;
  string str_spp_occurrence_filepath;
  string str_filename;
  string str_Env1Grids_filepath;
  string str_Env1_grids_namefilepath;
  string s_ouputfile_name;
  string str_species_number;
  string str_locations_filename;
  string str_locations_filepath;
  string str_spp_dispersal_filename;
  string str_condition_filename;
  string str_condition_filepath;
  string str_spp_parameter_filepath;
  string str_spp_occurrence_filename;
  string str_env_number;
  string str_spp_dispersal_filepath;
  string str_Vp_grid_filepath;
  string str_Keep_time;
  string str_hdr_name;
  string str_flt_name;

  // String arrays
  string str_env_grids_folder_path;  // To hold env folder paths: Max 10 environmental variables
  string str_env_grids_namefile;    // To hold env namefiles: Max 10 environmental variables
  string str_env_grids_files[10000];   // to hold the env grids filenames
  string str_spp_parameter_files[200];    // To hold species parameter files: Max 200 species
  string str_Z_grid_filename[50];    // To hold species Z namefiles: Max 100 environmental variables
  string str_Vp_grid_filename[50];    // To hold species Vp namefiles: Max 100 environmental variables

  // other hard-coded parameter arrays
  int b_Z_grid[50];
  int b_Vp_grid[50];
  int b_env_low_adaptation[50];
  int b_env_high_adaptation[50];
  float f_env_low_fund_limit[50];
  float f_env_high_fund_limit[50];
  double d_env_Z[50];
  double d_env_Vp[50];
  double d_spp_H2[50];
  double d_adapt_limit[50];
  double d_plasticity[50];
  double d_nat_sel_slope[50];

  // Declare pointers to local arrays
  short** s_condition;
  short** s_initial_occurrence;
  int** i_row_col;
  int** i_rel_row_col;
  double* d_revised_Z_upper_limit;
  double* d_revised_Vp_upper_limit;
  double* d_rel_prob;
  double* d_mean_Z_env;
  double* d_mean_Vp_env;

  // other pointers
  float* f_env_array;
  float* f_ptr;
  double* d_Z_env_array;
  double* d_Vp_env_array;

  TRaft* Raft; // pointer to the Raft

  TGridCell_ptr FromCell; // pointer to the From Cell (the focal cell)
  TGridCell_ptr ToCell; // pointer to the To Cell (the focal cell)

//// Open parameter file and read in the parameters/////////////////////////////
  ifstream ParamFile;
  ParamFile.open(arg_filename.c_str());
  if(ParamFile.is_open())
	{
	ParamFile >> str_buffer; // "NCOLS"
	ParamFile >> n_cols;
	ParamFile >> str_buffer; // "NROWS"
	ParamFile >> n_rows;
//OLD	ParamFile >> str_buffer; // "INPUT_FOLDER"
//OLD	i_buf=ParamFile.get(); // read blank space
//OLD	getline(ParamFile,str_buffer); // read the input folder
//OLD	str_input_folder_path=str_buffer.c_str(); //put the input folder name where it belongs
	ParamFile >> str_buffer; // "OUTPUT_FOLDER"
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read the folder path
	str_output_folder_path=str_buffer.c_str(); //put the folder name where it belongs
	ParamFile >> str_buffer; // "OUTPUT_ROOTNAME"
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read the name
	str_output_rootname=str_buffer.c_str(); //put the name where it belongs
	ParamFile >> str_buffer; // "WRITE_GRIDS_EACH_STEP"
	ParamFile >> b_write_each_step;
	ParamFile >> str_buffer; // "NUM_TIME_POINTS"
	ParamFile >> n_timepoints;
	ParamFile >> str_buffer; // "NUM_ENVIRONMENT_VARIABLES"
	ParamFile >> n_environments;
	ParamFile >> str_buffer; // "ENV_GRIDS_FOLDER"
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read rest of line
	str_env_grids_folder_path=str_buffer.c_str();
	ParamFile >> str_buffer; // "ENV_GRIDS_NAMEFILE"
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read rest of line
	str_env_grids_namefile=str_buffer.c_str();
//OLD	ParamFile >> str_buffer; // "SPECIES_GRIDS_FOLDER"
//OLD	i_buf=ParamFile.get(); // read blank space
//OLD	getline(ParamFile,str_buffer); // read the folder path
//OLD	str_spp_grids_folder_path=str_buffer.c_str(); //put the folder name where it belongs
	ParamFile >> str_buffer; // "SPECIES_OCCURRENCE_GRID"
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read the input filename
	str_spp_occurrence_filename=str_buffer.c_str(); //put the input filename where it belongs      //THIS IS NOW THE CONDITION AND THE INITIAL OCCURRENCE GRID
	ParamFile >> str_buffer; //MIN_SURVIVAL_PERCENTAGE
	ParamFile >> d_spp_survival_min;
	ParamFile >> str_buffer; //RESIDENT_POPULATION_ATTRIBUTE_WEIGHTING
	ParamFile >> d_resident_weighting;
	ParamFile >> str_buffer; //DISPERSAL_NEIGHBOURHOOD_FILE
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read the input folder
	str_spp_dispersal_filename=str_buffer.c_str(); //put the input folder name where it belongs
	ParamFile >> str_buffer; //SPECIES_LOCATION_FILE
	i_buf=ParamFile.get(); // read blank space
	getline(ParamFile,str_buffer); // read the input folder
	str_locations_filename=str_buffer.c_str(); //put the input folder name where it belongs
	for(i_env=0;i_env<n_environments;i_env++)
	  {
	  ParamFile >> str_buffer; //ENV_VARIABLE_X
	  ParamFile >> str_buffer; //LOWFUNDLIMIT
	  ParamFile >> f_env_low_fund_limit[i_env];
	  ParamFile >> str_buffer; //HIGHFUNDLIMIT
	  ParamFile >> f_env_high_fund_limit[i_env];
	  ParamFile >> str_buffer; //LOWADAPTATION
	  ParamFile >> b_env_low_adaptation[i_env];
	  ParamFile >> str_buffer; //HIGHADAPTATION
	  ParamFile >> b_env_high_adaptation[i_env];
	  if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		{
		ParamFile >> str_buffer; //ADAPTATIONLIMIT
		ParamFile >> d_adapt_limit[i_env];
		ParamFile >> str_buffer; //HERITABILITY
		ParamFile >> d_spp_H2[i_env];
		ParamFile >> str_buffer; //FITNESSCOST
		ParamFile >> d_nat_sel_slope[i_env];
		ParamFile >> str_buffer; //Z_GRID
		ParamFile >> b_Z_grid[i_env];  // 0=just use fund limit, 1=read in grid
		ParamFile >> str_buffer; //Z_GRIDNAME_OR_Z_VALUE
		if(b_Z_grid[i_env]>0)
		  {
		  i_buf=ParamFile.get(); // read blank space
		  getline(ParamFile,str_buffer); // read the input folder
		  str_Z_grid_filename[i_env]=str_buffer.c_str(); //put the input folder name where it belongs
		  } // end if b_Z_grid[i_env]>0
		else
		  ParamFile >> d_env_Z[i_env];
		ParamFile >> str_buffer; //VP_GRID
		ParamFile >> b_Vp_grid[i_env];  // 0=just use fund limit, 1=read in grid
		ParamFile >> str_buffer; //VP_GRIDNAME_OR_VP_VALUE
		if(b_Vp_grid[i_env]>0)
		  {
		  i_buf=ParamFile.get(); // read blank space
		  getline(ParamFile,str_buffer); // read the input folder
		  str_Vp_grid_filename[i_env]=str_buffer.c_str(); //put the input folder name where it belongs
		  } // end if b_Vp_grid[i_env]>0
		else
		  ParamFile >> d_env_Vp[i_env];
		ParamFile >> str_buffer; //PLASTICITY
		ParamFile >> d_plasticity[i_env];
		if(b_env_low_adaptation[i_env]>0)  // make corrections for plasticity
		  {
		  f_env_high_fund_limit[i_env] = f_env_high_fund_limit[i_env] + (float) (d_plasticity[i_env]);
		  d_plasticity[i_env] = d_plasticity[i_env] * -1;  // so we can just use one line to correct the environment for plasticity later on
		  } // end if b_env_low_adaptation[i_env]>0
		else
		  f_env_low_fund_limit[i_env] = f_env_low_fund_limit[i_env] - (float) (d_plasticity[i_env]);
		} //end if b_low_adaptation>0 || b_high_adaptation>0
	  } // end for i_env
	ParamFile >> str_buffer; // "PERCENTAGE_SUCCESS"
	ParamFile >> i_check;
	ParamFile.close();
	if(i_check != 100)  // ensure we have read in the parameters correctly. If not, get the hell out of there.
	  return -2;
	}	// end if opened
	else
	  return -1;
///////////////////////////////////////////////////////end parameter file read//

 //TEMP
  string str_checker_name = "C:/Users/mok010/Karels Files/Mokany Files/My R applications/AdaptR_package/AdaptR/data/outputs/AdaptR_CHECKER.txt";
  ofstream CHECKER(str_checker_name.c_str());
  CHECKER << "heritability " <<  d_spp_H2[0] << "\n";
  CHECKER << "fitness_cost " <<  d_nat_sel_slope[0] << "\n";
  CHECKER << "threshold " <<  d_adapt_limit[0] << "\n";
  CHECKER << "threshold_after_plasticity " << f_env_high_fund_limit[0]  << "\n";
  CHECKER << "plasticity " << d_plasticity[0]  << "\n";
  CHECKER << "low_adaptation " << b_env_low_adaptation[0]<< "\n";
  CHECKER << "high_adaptation " << b_env_high_adaptation[0]<< "\n";
  CHECKER.close();
// end TEMP  
  
  // set up difftime function to determine total processing time
  time_t start1, end1;
  start1 = time(NULL);

  // set the seed //////////////////////////////////////////////////////////////
  // use the parameter file name as a start
  str_buffer = arg_filename;
  i_number = 0;
  for(i_rel=0; i_rel< str_buffer.length(); i_rel++)
	{
	i_buf = int (str_buffer[i_rel]);
	i_number += i_buf;
	}
  // then add it to the time
  srand(time(NULL)*i_number);
  ///////////////////////////////////////////////////////////end set the seed //

  // Allocate arrays
  s_condition=AllocateShort2DArray(n_rows,
								   n_cols);
  s_initial_occurrence=AllocateShort2DArray(n_rows,
											n_cols);
  d_revised_Z_upper_limit = AllocateDouble1DArray(i_n_segments);
  d_revised_Vp_upper_limit = AllocateDouble1DArray(i_n_segments);
  d_mean_Z_env= AllocateDouble1DArray(n_environments);
  d_mean_Vp_env= AllocateDouble1DArray(n_environments);

  // open the condition grid and read it into a condition array/////////////////
//OLD//  str_spp_occurrence_filepath = str_spp_grids_folder_path + str_spp_occurrence_filename;
  str_spp_occurrence_filepath = str_spp_occurrence_filename; //NEW//
  ifstream Condition_FileStream;
  Condition_FileStream.open(str_spp_occurrence_filepath.c_str());
  // read in the header
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>n_totalcols;
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>n_totalrows;
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>d_xllcorner;
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>d_yllcorner;
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>d_cellsize;
  Condition_FileStream>>str_buffer;
  Condition_FileStream>>i_nodata;
  // Read in the data
  n_occupied = 0;
  for(i_row=0; i_row<n_rows; i_row++)
	{
	for(i_col=0; i_col<n_cols; i_col++)
	  {
	  Condition_FileStream >> f_data;
	  if(f_data > -9000)
		{
		s_condition[i_row][i_col] = 1;
		if(f_data > 0)
		  {
		  s_initial_occurrence[i_row][i_col] = 1;
		  n_occupied += 1;
		  }
		else
		  s_initial_occurrence[i_row][i_col] = (-9999);
		} // end if f_data>-9000
	  else
		{
		s_condition[i_row][i_col] = (-9999);
		s_initial_occurrence[i_row][i_col] = (-9999);
		} // end else
	  } // end for i_col
	} // end for i_row
  Condition_FileStream.close();
////////////////////////////////////////////////////end read in condition grid//

  // make a new raft ///////////////////////////////////////////////////////////
  Raft=new TRaft(n_rows,
				 n_cols,
				 s_condition);

  // Allocate the right amount of space on the raft
  Raft->LoadEmptyRaftData(n_rows,
						  n_cols,
						  n_environments);
/////////////////////////////////////////////////////////////end make new raft//

  // Fill the relative arrays to allow acounting for when the distribution of
  // attribute values crosses the upper critical threshold
  i_index_for_mean = (i_n_segments - 1) / 2;
  d_increment_per_index = i_index_for_mean / d_variance_limit;
  Calculate_upper_threshold_correction(i_n_segments,
									   i_index_for_mean,
									   d_variance_limit,
									   d_revised_Z_upper_limit,		// gets filled here
									   d_revised_Vp_upper_limit);   // gets filled here

  // Read in the environment files to the filename array ///////////////////////
	// open the name file
	str_buffer=str_env_grids_namefile;
	str_filename=str_buffer.c_str();
//OLD//	str_Env1_grids_namefilepath = str_input_folder_path + str_filename;
	str_Env1_grids_namefilepath = str_filename; //NEW//
	ifstream Env1Grids_NameFileStream;
	Env1Grids_NameFileStream.open(str_Env1_grids_namefilepath.c_str());
	// read the env grids folder
	str_buffer=str_env_grids_folder_path;
	str_env1_grids_folder_path = str_buffer.c_str();
	for(i_time=0; i_time<n_timepoints; i_time++)
	  {
	  getline(Env1Grids_NameFileStream,str_buffer); // read the input folder
	  str_filename=str_buffer.c_str(); //put the input folder name where it belongs
	  str_Env1_grids_namefilepath = str_env1_grids_folder_path + str_filename + ".mcg";
	  str_env_grids_files[i_time] = str_Env1_grids_namefilepath;
	  } // end for i_time
	Env1Grids_NameFileStream.close();
  ///////////////////////////////////////////////////end read in Env filenames//


	// Open the species occurrence grid and read it into the raft //////////////
//OLD//	str_spp_occurrence_filepath = str_spp_grids_folder_path + str_spp_occurrence_filename;
    str_spp_occurrence_filepath = str_spp_occurrence_filename; //NEW//
	ifstream SppOccurrence_FileStream;
	SppOccurrence_FileStream.open(str_spp_occurrence_filepath.c_str());
	// read in the header
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>n_totalcols;
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>n_totalrows;
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>d_xllcorner;
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>d_yllcorner;
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>d_cellsize;
	SppOccurrence_FileStream>>str_buffer;
	SppOccurrence_FileStream>>i_nodata;
	for(i_row=0; i_row<n_rows; i_row++)
	  {
	  for(i_col=0; i_col<n_cols; i_col++)
		{
		SppOccurrence_FileStream >> f_data;
		FromCell = Raft->GetCell(i_row,
								 i_col);
		if(FromCell != NULL) // Check if this cell is land
		  {
		  if(f_data>0)
			{
			FromCell->SetOccurrence(s_one);
			FromCell->SetTimeLastPresent(s_zero);
			FromCell->SetNewlyColonised(s_zero);
			} // end if f_data>0
		  else
			{
			FromCell->SetOccurrence(s_zero);
			FromCell->SetTimeLastPresent(s_nodata);
			FromCell->SetNewlyColonised(s_zero);
			} // end else f_data>0
		  } // end if FromCell != NULL
		} // end for i_col
	  } // end for i_row
	SppOccurrence_FileStream.close();
	////////////////////////////////////////////////////end read in occurrence//

	// Read in the dispersal neighbourhood file (.dna) /////////////////////////
//OLD//	str_spp_dispersal_filepath = str_spp_grids_folder_path + str_spp_dispersal_filename;
    str_spp_dispersal_filepath = str_spp_dispersal_filename;//NEW//
	ifstream Spp_Disp_FileStream;
	Spp_Disp_FileStream.open(str_spp_dispersal_filepath.c_str());
	Spp_Disp_FileStream>>n_rel;
	// allocate the arrays
	i_rel_row_col=AllocateInt2DArray(n_rel,
									 2);
	d_rel_prob=AllocateDouble1DArray(n_rel);
	for(i_rel=0; i_rel<n_rel; i_rel++)
	  {
	  Spp_Disp_FileStream>>i_rel_row_col[i_rel][0];
	  Spp_Disp_FileStream>>i_rel_row_col[i_rel][1];
	  Spp_Disp_FileStream>>d_rel_prob[i_rel];
	  } // end for i_rel
	Spp_Disp_FileStream.close();
	////////////////////////////////////////////////end read in dispersal file//
	
	// Read in the locations to report on///////////////////////////////////////
	// Open the file holding the locations to report on
//OLD//	str_locations_filepath = str_spp_grids_folder_path + str_locations_filename;
    str_locations_filepath = str_locations_filename;//NEW//	
	ifstream SppLocations_FileStream;
	SppLocations_FileStream.open(str_locations_filepath.c_str());
	// read in the header
	SppLocations_FileStream >> str_buffer; //"n_locations"
	SppLocations_FileStream >> n_locations;
	SppLocations_FileStream >> str_buffer; //"Longitude"
	SppLocations_FileStream >> str_buffer; //"Latitude"
	// Allocate and array for the row & col
	i_row_col=AllocateInt2DArray(n_locations,
								 2);
	// Loop through the locations, write the lat & long to output file, calculate
	// the row and column and put it in the array
	d_yulcorner = d_yllcorner + (double) (n_rows * d_cellsize);
	for(i_loc=0; i_loc<n_locations; i_loc++)
	  {
	  SppLocations_FileStream >> d_pt_Longitude; //"Longitude"
	  SppLocations_FileStream >> d_pt_Latitude; //"Latitude"
	  i_row_col[i_loc][0] = (int) ((d_yulcorner - d_pt_Latitude) / d_cellsize);  // row
	  i_row_col[i_loc][1] = (int) ((d_pt_Longitude - d_xllcorner) / d_cellsize); // col
	  } // end for i_loc
	SppLocations_FileStream.close();
	// Open the output file to write the location results to
	s_ouputfile_name = str_output_folder_path+str_output_rootname+"_species_locations_results.txt";
	ofstream LocationsOutputFileStream(s_ouputfile_name.c_str());
	LocationsOutputFileStream << "Time Location Species_Present ";
	for(i_env=1;i_env<(n_environments+1);i_env++)
	  LocationsOutputFileStream << "Value_Env_" << i_env << " ";
	for(i_env=1;i_env<(n_environments+1);i_env++)
	  LocationsOutputFileStream << "Z_Env_" << i_env << " ";
	for(i_env=1;i_env<(n_environments+1);i_env++)
	  LocationsOutputFileStream << "Vp_Env_" << i_env << " ";
	LocationsOutputFileStream << endl;
	////////////////////////////////////////////////end read in locations file//

	// Initialise the output file to report on the area occupied at each time point
	s_ouputfile_name = str_output_folder_path+str_output_rootname+"_occupancy_summary.txt";
	ofstream SumOutputFileStream(s_ouputfile_name.c_str());
	SumOutputFileStream << "Time Cells_Occupied Initial_Cells_Occupied Not_initial_Cells_Occupied ";
	for(i_env=1;i_env<(n_environments+1);i_env++)
	  {
	  if(b_env_low_adaptation[(i_env-1)]>0 || b_env_high_adaptation[(i_env-1)]>0)
		SumOutputFileStream << "Mean_Z_Env_" << i_env << " " << "Mean_Vp_Env_" << i_env << " ";
	  } // end for i_env
	SumOutputFileStream << endl;
	SumOutputFileStream << "0" << " " << n_occupied << " " << n_occupied << " " << "0" << " ";
	/////////////////////////////////// end initialise area of occupancy file //


	// Initialise relevent Z & Vp grids ////////////////////////////////////////
	for(i_env=0;i_env<n_environments;i_env++)
	  {
	  if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		{
		if(b_Z_grid[i_env]>0)
		  {
		  // read in the Z grid
		  str_buffer=str_Z_grid_filename[i_env];
		  str_filename=str_buffer.c_str();
//OLD//		  str_spp_occurrence_filepath = str_spp_grids_folder_path + str_filename;
          str_spp_occurrence_filepath = str_filename;
		  ifstream Spp_Z_FileStream;
		  Spp_Z_FileStream.open(str_spp_occurrence_filepath.c_str());
		  // read in the header
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>n_totalcols;
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>n_totalrows;
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>d_xllcorner;
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>d_yllcorner;
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>d_cellsize;
		  Spp_Z_FileStream>>str_buffer;
		  Spp_Z_FileStream>>i_nodata;
		  d_mean_Z_env[0] = 0;
		  for(i_row=0; i_row<n_rows; i_row++)
			{
			for(i_col=0; i_col<n_cols; i_col++)
			  {
			  Spp_Z_FileStream >> d_data;
			  FromCell = Raft->GetCell(i_row,
									   i_col);
			  if(FromCell != NULL) // Check if this cell is land
				{
				s_present = FromCell->GetOccurrence();
				d_Z_env_array = FromCell->GetCurrentZ();
				if(s_present>0)
				  {
				  d_Z_env_array[i_env] = d_data;
				  d_mean_Z_env[0] += d_data;
				  }
				else
				  d_Z_env_array[i_env] = d_nodata;
				} // end if FromCell != NULL
			  } // end for i_col
			} // end for i_row
		  Spp_Z_FileStream.close();
		  if(n_occupied>0)
			d_mean_Z_env[0] = d_mean_Z_env[0] / (double) (n_occupied);
		  SumOutputFileStream << d_mean_Z_env[0] << " ";
		  } // end if b_Z_grid[i_env]>0
		else
		  {
		  // set all the grid values to the specified fundamental limit value
		  for(i_row=0; i_row<n_rows; i_row++)
			{
			for(i_col=0; i_col<n_cols; i_col++)
			  {
			  FromCell = Raft->GetCell(i_row,
									   i_col);
			  if(FromCell != NULL) // Check if this cell is land
				{
				s_present = FromCell->GetOccurrence();
				d_Z_env_array = FromCell->GetCurrentZ();
				if(s_present>0)
				  d_Z_env_array[i_env] = d_env_Z[i_env];
				else
				  d_Z_env_array[i_env] = d_nodata;
				} // end if FromCell != NULL
			  } // end for i_col
			} // end for i_row
		  SumOutputFileStream << d_env_Z[i_env] << " ";
		  } // end else b_Z_grid[i_env]>0
		if(b_Z_grid[i_env]>0)
		  {
		  // read in the Vp grid
		  str_buffer=str_Vp_grid_filename[i_env];
		  str_filename=str_buffer.c_str();
//OLD//		  str_Vp_grid_filepath = str_spp_grids_folder_path + str_filename;
          str_Vp_grid_filepath = str_filename;//NEW//
		  ifstream Spp_Vp_FileStream;
		  Spp_Vp_FileStream.open(str_Vp_grid_filepath.c_str());
		  // read in the header
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>n_totalcols;
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>n_totalrows;
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>d_xllcorner;
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>d_yllcorner;
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>d_cellsize;
		  Spp_Vp_FileStream>>str_buffer;
		  Spp_Vp_FileStream>>i_nodata;
		  d_mean_Vp_env[0] = 0;
		  for(i_row=0; i_row<n_rows; i_row++)
			{
			for(i_col=0; i_col<n_cols; i_col++)
			  {
			  Spp_Vp_FileStream >> d_data;
			  FromCell = Raft->GetCell(i_row,
									   i_col);
			  if(FromCell != NULL) // Check if this cell is land
				{
				s_present = FromCell->GetOccurrence();
				d_Vp_env_array = FromCell->GetCurrentVp();
				if(s_present>0)
				  {
				  d_Vp_env_array[i_env] = d_data;
				  d_mean_Vp_env[0] += d_data;
				  } // end if s_present>0
				else
				  d_Vp_env_array[i_env] = d_nodata;
				} // end if FromCell != NULL
			  } // end for i_col
			} // end for i_row
		  Spp_Vp_FileStream.close();
		  if(n_occupied>0)
			d_mean_Vp_env[0] = d_mean_Vp_env[0] / (double) (n_occupied);
		  SumOutputFileStream << d_mean_Vp_env[0] << " ";
		  } // end if b_Z_grid[i_env]>0
		else
		  {
		  // set all the grid values to the specified fundamental limit value
		  for(i_row=0; i_row<n_rows; i_row++)
			{
			for(i_col=0; i_col<n_cols; i_col++)
			  {
			  FromCell = Raft->GetCell(i_row,
									   i_col);
			  if(FromCell != NULL) // Check if this cell is land
				{
				s_present = FromCell->GetOccurrence();
				d_Vp_env_array = FromCell->GetCurrentVp();
				if(s_present>0)
				  d_Vp_env_array[i_env] = d_env_Vp[i_env];
				else
				  d_Vp_env_array[i_env] = d_nodata;
				} // end if FromCell != NULL
			  } // end for i_col
			} // end for i_row
		  SumOutputFileStream << d_env_Vp[i_env] << " ";
		  } // end else b_Vp_grid[i_env]>0
		} // end if (b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
	  } // end for i_env
      SumOutputFileStream << endl;
	  //////////////////////////////////////////// end initialise Z & Vp data //

	// Now loop through each time point //~~~~////~~~~////~~~~////~~~~////~~~~//
	for(i_time=0; i_time<n_timepoints; i_time++)
	  {
	  // Read in the environment data for i_time ///////////////////////////////
	  str_buffer=str_env_grids_files[i_time];
	  str_Env1Grids_filepath=str_buffer.c_str();
	  i_data_read = Raft->ReadInWholeENVStack(str_Env1Grids_filepath,
											  n_environments);
	  ////////////////////////////end read in the environment data for i_time //

	  // Write the current Z to the locations output file //////////////////////
	  for(i_loc=0; i_loc<n_locations; i_loc++)
		{
		i_number = i_loc + 1;
		i_row = i_row_col[i_loc][0];
		i_col = i_row_col[i_loc][1];
		FromCell = Raft->GetCell(i_row,
								 i_col);
		if(FromCell != NULL) // Check if this cell is land
		  {
		  LocationsOutputFileStream << i_time << " ";
		  LocationsOutputFileStream << i_number << " ";
		  s_occurrence = FromCell->GetOccurrence();
		  LocationsOutputFileStream << s_occurrence << " ";
		  f_env_array = FromCell->GetEnvData();
		  for(i_env=0;i_env<n_environments;i_env++)
			LocationsOutputFileStream << f_env_array[i_env] << " ";
		  d_Z_env_array = FromCell->GetCurrentZ();
		  for(i_env=0;i_env<n_environments;i_env++)
			LocationsOutputFileStream << d_Z_env_array[i_env] << " ";
		  d_Vp_env_array = FromCell->GetCurrentVp();
		  for(i_env=0;i_env<n_environments;i_env++)
			LocationsOutputFileStream << d_Vp_env_array[i_env] << " ";
		  LocationsOutputFileStream << endl;
		  } // end if FromCell != NULL
		} // end for i_loc
      //////////////////////////////////////////////end write location outputs//

	  // Add plasticity offset to the relevant environment layer ///////////////
	  for(i_env=0;i_env<n_environments;i_env++)
		{
		if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		  {
		  if(d_plasticity[i_env]>0 || d_plasticity[i_env]<0)
			{
            for(i_row=0; i_row<n_rows; i_row++)
			  {
			  for(i_col=0; i_col<n_cols; i_col++)
				{
				FromCell = Raft->GetCell(i_row,
										 i_col);
				if(FromCell != NULL) // Check if this cell is land
				  {
				  f_env_array = FromCell->GetEnvData();
				  f_env_array[i_env] = f_env_array[i_env] - (float) (d_plasticity[i_env]);  // we have already adjusted adaptation at low limit so -= will work either way
				  } // end if FromCell != NULL
				} // end for i_col
			  } // end for i_row
			} // end if d_plasticity[i_env]>0
		  } // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
		} // end for i_env
	  ///////////////////////////////////////////////end add plasticity offset//

	  // Now loop through the grid to determine persistence of occupied grid cells//
	  Assess_Persistence_and_Adaptation(Raft,
										n_rows,
										n_cols,
										n_environments,
										i_time,
										i_index_for_mean,
										i_n_segments,
										d_spp_survival_min,
										d_increment_per_index,
										b_env_low_adaptation,
										b_env_high_adaptation,
										f_env_low_fund_limit,
										f_env_high_fund_limit,
										d_adapt_limit,
										d_revised_Z_upper_limit,
										d_revised_Vp_upper_limit,
										d_spp_H2,
										d_nat_sel_slope);
	  /////////////////////////////////////end assess persistence & adaptation//

	  // Now see if there's any dispersal to unoccupied grid cells /////////////
	  Implement_Dispersal_To_Unoccupied_Cells(Raft,
											  n_rows,
											  n_cols,
											  n_rel,
											  n_environments,
											  b_env_low_adaptation,
											  b_env_high_adaptation,
											  i_rel_row_col,
											  f_env_low_fund_limit,
											  f_env_high_fund_limit,
											  d_rel_prob,
											  d_adapt_limit);
	  ///////////////////////////////////////end dispersal to uniccupied cells//

	  // Now account for change in adaptive parameters due to dispersal amongst
	  // occupied grid cells
	  Implement_Adaptive_Parameter_Averaging_Mixture(Raft,
													 n_rows,
													 n_cols,
													 n_rel,
													 n_environments,
													 d_resident_weighting,
													 b_env_low_adaptation,
													 b_env_high_adaptation,
													 i_rel_row_col,
													 d_rel_prob,
													 d_adapt_limit);
	  ////////////////////////////////////// end adaptive parameter averaging //

	  // Now make all newly colonised grid cells actually occupied//////////////
	  convert_colonised_cellls_to_occupied(Raft,
										   n_rows,
										   n_cols);
	  ////////////////////////////make all newly colonised grid cells occupied//

	  // Now calculate the occurrence & Mean Z & Vp summary & write it to file //
	  Generate_occurrence_summary(Raft,
								  n_rows,
								  n_cols,
								  n_environments,
								  b_env_low_adaptation,
								  b_env_high_adaptation,
								  s_initial_occurrence,
								  d_mean_Z_env,
								  d_mean_Vp_env,
								  &n_occupied,
								  &n_initial_occupied);

	  i_now = i_time + 1;
	  n_not_initial_occupied = n_occupied - n_initial_occupied;
	  SumOutputFileStream << i_now << " " << n_occupied << " " << n_initial_occupied << " " << n_not_initial_occupied << " ";
	  for(i_env=0;i_env<n_environments;i_env++)
		{
		if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		  SumOutputFileStream << d_mean_Z_env[i_env] << " " << d_mean_Vp_env[i_env] << " ";
		} // end for i_env
	  SumOutputFileStream << endl;
	  ////////////////////////////////////////////// End calc occurence summary//

	  // NEW - 18 Dec 2015 //
	  // Write time T occurrence and Z to float grid ///////////////////////////
	  if(b_write_each_step > 0)
		{
	  stringstream TimeStrStream;
	  TimeStrStream << i_time;
	  str_Keep_time = TimeStrStream.str();

	  ofstream HeaderFileStream;
      ofstream FltFileStream;

	  // Occurrence
	  str_hdr_name = str_output_folder_path+str_output_rootname+"_occurrence_T"+str_Keep_time+".hdr";
	  str_flt_name = str_output_folder_path+str_output_rootname+"_occurrence_T"+str_Keep_time+".flt";
	  HeaderFileStream.open(str_hdr_name.c_str());
	  FltFileStream.open(str_flt_name.c_str(), ios::out |ios::binary);

	  HeaderFileStream<<"ncols"<<"\t\t"<<n_cols<<endl;
	  HeaderFileStream<<"nrows"<<"\t\t"<<n_rows<<endl;
	  HeaderFileStream<<"xllcorner"<<"\t"<<d_xllcorner<<endl;
	  HeaderFileStream<<"yllcorner"<<"\t"<<d_yllcorner<<endl;
	  HeaderFileStream<<"cellsize"<<"\t"<<d_cellsize<<endl;
	  HeaderFileStream<<"NODATA_value"<<"\t"<<i_nodata<<endl;
	  HeaderFileStream<<"BYTEORDER\tLSBFIRST"<<endl;
	  HeaderFileStream.close();

	  for(i_row=0; i_row<n_rows; i_row++)
		{
		for(i_col=0; i_col<n_cols; i_col++)
		  {
		  FromCell = Raft->GetCell(i_row,
								 i_col);
		  if(FromCell != NULL) // Check if this cell is land
			{
			s_present_time = FromCell->GetOccurrence();
			if(s_present_time>0)
			  {
			  f_ptr = &(f_one);
			  FltFileStream.write((char *)f_ptr, sizeof(float));
			  } // end if s_present_time>0
			else
			  {
			  f_ptr = &(f_nodata);
			  FltFileStream.write((char *)f_ptr, sizeof(float));
			  } // end else
			}
		  else
			{
			f_ptr = &(f_nodata);
			FltFileStream.write((char *)f_ptr, sizeof(float));
			}
		  } // end for i_col
		} // end for i_row
	  FltFileStream.close();

	  // Z
	  for(i_env=0;i_env<n_environments;i_env++)
		{
		if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		  {
		  i_this_env = i_env+1;
		  stringstream EnvStrStream;
		  EnvStrStream << i_this_env;
		  str_env_number = EnvStrStream.str();
		  // write the reults to file - final z value ////////////////////////////
		  str_hdr_name = str_output_folder_path+str_output_rootname+"_Env_"+str_env_number+"_Z_T"+str_Keep_time+".hdr";
		  str_flt_name = str_output_folder_path+str_output_rootname+"_Env_"+str_env_number+"_Z_T"+str_Keep_time+".flt";
		  HeaderFileStream.open(str_hdr_name.c_str());
		  FltFileStream.open(str_flt_name.c_str(), ios::out |ios::binary);

		  HeaderFileStream<<"ncols"<<"\t\t"<<n_cols<<endl;
		  HeaderFileStream<<"nrows"<<"\t\t"<<n_rows<<endl;
		  HeaderFileStream<<"xllcorner"<<"\t"<<d_xllcorner<<endl;
		  HeaderFileStream<<"yllcorner"<<"\t"<<d_yllcorner<<endl;
		  HeaderFileStream<<"cellsize"<<"\t"<<d_cellsize<<endl;
		  HeaderFileStream<<"NODATA_value"<<"\t"<<i_nodata<<endl;
		  HeaderFileStream<<"BYTEORDER\tLSBFIRST"<<endl;
		  HeaderFileStream.close();

		for(i_row=0; i_row<n_rows; i_row++)
		  {
		  for(i_col=0; i_col<n_cols; i_col++)
			{
			FromCell = Raft->GetCell(i_row,
									 i_col);
			if(FromCell != NULL) // Check if this cell is land
			  {
			  s_present_time = FromCell->GetOccurrence();
			  if(s_present_time>0)
				{
				d_Z_env_array = FromCell->GetCurrentZ();
				f_Z = float(d_Z_env_array[i_env]);
				f_ptr = &(f_Z);
				FltFileStream.write((char *)f_ptr, sizeof(float));
				} // end if s_present_time>0
			  else
				{
				f_ptr = &(f_nodata);
				FltFileStream.write((char *)f_ptr, sizeof(float));
				}
			  } // end if FromCell != NULL
			else
			  {
			  f_ptr = &(f_nodata);
			  FltFileStream.write((char *)f_ptr, sizeof(float));
			  }
			} // end for i_col
		  } // end for i_row
		FltFileStream.close();
		} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
	  } // end for i_env
	  } // end if write each step
	 //////////////////////// end Write time T occurrence and Z to float grid //


	  //////////////////////////////////////////////////////////////////////////

	  } // end for i_time

	// Write the reults to file - final occurrence /////////////////////////////
	s_ouputfile_name = str_output_folder_path+str_output_rootname+"_final_occurrence.asc";
	ofstream OutputFileStream(s_ouputfile_name.c_str());
	OutputFileStream<<"ncols         ";
	OutputFileStream<<n_cols<<endl;
	OutputFileStream<<"nrows         ";
	OutputFileStream<<n_rows<<endl;
	OutputFileStream<<"xllcorner     ";
	OutputFileStream<<d_xllcorner<<endl;
	OutputFileStream<<"yllcorner     ";
	OutputFileStream<<d_yllcorner<<endl;
	OutputFileStream<<"cellsize      ";
	OutputFileStream<<d_cellsize<<endl;
	OutputFileStream<<"NODATA_value  ";
	OutputFileStream<<i_nodata<<endl;
	for(i_row=0; i_row<n_rows; i_row++)
	  {
	  for(i_col=0; i_col<n_cols; i_col++)
		{
		FromCell = Raft->GetCell(i_row,
								 i_col);
		if(FromCell != NULL) // Check if this cell is land
		  {
		  s_present_time = FromCell->GetOccurrence();
		  OutputFileStream << s_present_time << " ";
		  }
		else
		  OutputFileStream << i_nodata << " ";
		} // end for i_col
	  OutputFileStream << endl;
	  } // end for i_row
	OutputFileStream.close();
	////////////////////////////////////////end write final occurrence to file//

	for(i_env=0;i_env<n_environments;i_env++)
	  {
	  if(b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0)
		{
		i_this_env = i_env+1;
		stringstream EnvStrStream;
		EnvStrStream << i_this_env;
		str_env_number = EnvStrStream.str();
		// write the reults to file - final z value ////////////////////////////
		s_ouputfile_name = str_output_folder_path+str_output_rootname+"_Env_"+str_env_number+"_final_Z.asc";
		ofstream OutputFileStream2(s_ouputfile_name.c_str());
		OutputFileStream2<<"ncols         ";
		OutputFileStream2<<n_cols<<endl;
		OutputFileStream2<<"nrows         ";
		OutputFileStream2<<n_rows<<endl;
		OutputFileStream2<<"xllcorner     ";
		OutputFileStream2<<d_xllcorner<<endl;
		OutputFileStream2<<"yllcorner     ";
		OutputFileStream2<<d_yllcorner<<endl;
		OutputFileStream2<<"cellsize      ";
		OutputFileStream2<<d_cellsize<<endl;
		OutputFileStream2<<"NODATA_value  ";
		OutputFileStream2<<i_nodata<<endl;
		for(i_row=0; i_row<n_rows; i_row++)
		  {
		  for(i_col=0; i_col<n_cols; i_col++)
			{
			FromCell = Raft->GetCell(i_row,
									 i_col);
			if(FromCell != NULL) // Check if this cell is land
			  {
			  s_present_time = FromCell->GetOccurrence();
			  if(s_present_time>0)
				{
				d_Z_env_array = FromCell->GetCurrentZ();
				OutputFileStream2 << d_Z_env_array[i_env] << " ";
				} // end if s_present_time>0
			  else
				OutputFileStream2 << i_nodata << " ";
			  } // end if FromCell != NULL
			else
			  OutputFileStream2 << i_nodata << " ";
			} // end for i_col
		  OutputFileStream2 << endl;
		  } // end for i_row
		OutputFileStream2.close();
		///////////////////////////////////////end write final z-value to file//

		// write the reults to file - final Vp value ///////////////////////////
		s_ouputfile_name = str_output_folder_path+str_output_rootname+"_Env_"+str_env_number+"_final_Vp.asc";
		ofstream OutputFileStream3(s_ouputfile_name.c_str());
		OutputFileStream3<<"ncols         ";
		OutputFileStream3<<n_cols<<endl;
		OutputFileStream3<<"nrows         ";
		OutputFileStream3<<n_rows<<endl;
		OutputFileStream3<<"xllcorner     ";
		OutputFileStream3<<d_xllcorner<<endl;
		OutputFileStream3<<"yllcorner     ";
		OutputFileStream3<<d_yllcorner<<endl;
		OutputFileStream3<<"cellsize      ";
		OutputFileStream3<<d_cellsize<<endl;
		OutputFileStream3<<"NODATA_value  ";
		OutputFileStream3<<i_nodata<<endl;
		for(i_row=0; i_row<n_rows; i_row++)
		  {
		  for(i_col=0; i_col<n_cols; i_col++)
			{
			FromCell = Raft->GetCell(i_row,
									 i_col);
			if(FromCell != NULL) // Check if this cell is land
			  {
			  s_present_time = FromCell->GetOccurrence();
			  if(s_present_time>0)
				{
				d_Vp_env_array = FromCell->GetCurrentVp();
				OutputFileStream3 << d_Vp_env_array[i_env] << " ";
				} // end if s_present_time>0
			  else
                OutputFileStream3 << i_nodata << " ";
			  } // end if FromCell != NULL
			else
			  OutputFileStream3 << i_nodata << " ";
			} // end for i_col
		  OutputFileStream3 << endl;
		  } // end for i_row
		OutputFileStream3.close();
		//////////////////////////////////////end write final Vp-value to file//
		} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
	  } // end for i_env

  // catch total processing time
  end1 = time(NULL);
  f_total_time = difftime(end1, start1);
  i_runtime = (int) (f_total_time);

  // Close the files
  LocationsOutputFileStream.close();
  SumOutputFileStream.close();

  // Free the arrays
  FreeShort2DArray(s_condition);
  FreeShort2DArray(s_initial_occurrence);
  FreeInt2DArray(i_rel_row_col);
  FreeDouble1DArray(d_rel_prob);
  FreeInt2DArray(i_row_col);
  FreeDouble1DArray(d_revised_Z_upper_limit);
  FreeDouble1DArray(d_revised_Vp_upper_limit);
  FreeDouble1DArray(d_mean_Z_env);
  FreeDouble1DArray(d_mean_Vp_env);

  // Free the allocated space on the raft
  Raft->FreeRaftData(n_rows,
					 n_cols);

  return(i_runtime);

} // end Run_Adaptor

//------------------------------------------------------------------------------
// A function that calculates an approximation to the standard normal cumulative
// distribution function.
// 'Borrowed' from : http://www.johndcook.com/cpp_phi.html
double phi(double arg_d_x)
{
  // Declare local variables
  int sign;
  double x;
  double t;
  double y;
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;


  // Save the sign of x
  sign = 1;
  if(arg_d_x < 0)
	sign = -1;
  x = fabs(arg_d_x)/sqrt(2.0);

  // A&S formula 7.1.26
  t = 1.0/(1.0 + p*x);
  y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

  return 0.5*(1.0 + sign*y);

} // end phi

//------------------------------------------------------------------------------
// A function that calculates survival percentage (s) given three attribute:
// the mean value of the normal distribution, the standard deviation, and the
// observed x-value.
double calculate_survival_percentage(double arg_d_x,
									 double arg_d_Z,
									 double arg_d_Vp)
{
  double d_relative_x;
  double d_relative_standard_x;
  double d_cummulative;
  double d_survival;

  d_relative_x = arg_d_x - arg_d_Z;
  d_relative_standard_x = d_relative_x / arg_d_Vp;
  d_cummulative = phi(d_relative_standard_x);
  d_survival = (1 - d_cummulative) * 100; // not sure if output is needed in fraction, or percent (x100)
  return(d_survival);

} // end calculate_survival_percentage

//------------------------------------------------------------------------------
// A function that calculates survival percentage (s) given three attribute:
// the mean value of the normal distribution, the standard deviation, and the
// observed x-value.
double calculate_survival_percentage_LowLimit(double arg_d_x,
											  double arg_d_Z,
											  double arg_d_Vp)
{
  double d_relative_x;
  double d_relative_standard_x;
  double d_cummulative;
  double d_survival;

  d_relative_x = arg_d_x - arg_d_Z;
  d_relative_standard_x = d_relative_x / arg_d_Vp;
  d_cummulative = phi(d_relative_standard_x);
  d_survival = d_cummulative * 100; // not sure if output is needed in fraction, or percent (x100)
  return(d_survival);

} // end calculate_survival_percentage_LowLimit

//------------------------------------------------------------------------------
// A function that calculates proportion of values above x given three attribute:
// the mean value of the normal distribution, the variance (converted to standard
// deviation within the function), and the observed x-value.
double calculate_above_limit_proportion(double arg_d_x,     //x
										double arg_d_Z,     //mean
										double arg_d_Vp)    //variance
{
  double d_relative_x;
  double d_relative_standard_x;
  double d_cummulative;
  double d_above;

  d_relative_x = arg_d_x - arg_d_Z;
  d_relative_standard_x = d_relative_x / sqrt(arg_d_Vp);
  d_cummulative = phi(d_relative_standard_x);
  d_above = (1 - d_cummulative); // not sure if output is needed in fraction, or percent (x100)
  return(d_above);

} // end calculate_above_limit_proportion
//------------------------------------------------------------------------------
// A function that calculates proportion of values below x given three attribute:
// the mean value of the normal distribution, the variance (converted to standard
// deviation within the function), and the observed x-value.
double calculate_below_limit_proportion(double arg_d_x,     //x
										double arg_d_Z,     //mean
										double arg_d_Vp)    //variance
{
  double d_relative_x;
  double d_relative_standard_x;
  double d_cummulative;

  d_relative_x = arg_d_x - arg_d_Z;
  d_relative_standard_x = d_relative_x / sqrt(arg_d_Vp);
  d_cummulative = phi(d_relative_standard_x);
  return(d_cummulative);

} // end calculate_below_limit_proportion

//------------------------------------------------------------------------------
// A function to calculate the relevent correction of the mean (Z) and variance (Vp)
// for when the upper critical threshold occurs within the specified bounds of the
// normal distribution of attribute values.
void Calculate_upper_threshold_correction(int arg_i_n_segments,
										  int arg_i_index_for_mean,
										  double arg_d_variance_limit,
										  double* arg_d_revised_Z_upper_limit,		// gets filled here
										  double* arg_d_revised_Vp_upper_limit)   // gets filled here
{
  // declare local variables
  int i_seg, i_incr;
  double d_potential_population = 100000;
  double d_actual_population;
  double d_population_against_threshold;
  double d_increment_per_index;
  double d_x, d_x_prior;
  double d_increment_proportion;
  double d_individuals;
  double d_mean_sum;
  double d_var_sum;

  // pointers for local arrays
  double* d_standardised_x;
  int* i_n_individuals;

  // allocate local arrays
  d_standardised_x = AllocateDouble1DArray(arg_i_n_segments);
  i_n_individuals = AllocateInt1DArray(arg_i_n_segments);

// TEMP - Write the filenames to a temporary file
/*
string str_output_filename;
str_output_filename = "C:\\temp\\Test_outputs\\Upper_correction.txt";
ofstream TempOutFileStream(str_output_filename.c_str());
*/
// END TEMP

  // calculate the increment per array index
  d_increment_per_index = arg_i_index_for_mean / arg_d_variance_limit;

  // Loop through and determine the standardised x-value for each segment of the
  // distribution
  for(i_seg=0; i_seg<arg_i_n_segments; i_seg++)
	d_standardised_x[i_seg] = (i_seg - arg_i_index_for_mean) / d_increment_per_index;

  // Now calculate the proportion of the distribution that sits within each segment
  // and also the number of individuals within each segment
  d_actual_population = 0;
  d_x = d_standardised_x[0];
  d_x += (1/(2*d_increment_per_index));
  d_increment_proportion = phi(d_x);
  d_individuals = d_increment_proportion * d_potential_population;
  i_n_individuals[0] = (int) (round_double(d_individuals));
  d_actual_population += i_n_individuals[0];
  for(i_seg=1; i_seg<arg_i_n_segments; i_seg++)
	{
	d_x = d_standardised_x[i_seg];
	d_x += (1/(2*d_increment_per_index));
	d_x_prior = d_standardised_x[(i_seg-1)];
	d_x_prior += (1/(2*d_increment_per_index));
	d_increment_proportion = phi(d_x) - phi(d_x_prior);
	d_individuals = d_increment_proportion * d_potential_population;
	i_n_individuals[i_seg] = (int) (round_double(d_individuals));
	d_actual_population += i_n_individuals[i_seg];
	} // end for i_seg

  // Now go through and calculate the mean and variance for when the critical
  // threshold is at each level through the distribution
  d_population_against_threshold = d_actual_population;
  arg_d_revised_Z_upper_limit[0] = d_standardised_x[0];
  arg_d_revised_Vp_upper_limit[0] = 0.00001; // something very small
  for(i_seg=1; i_seg<(arg_i_n_segments-1); i_seg++)
	{
	d_population_against_threshold = d_population_against_threshold - i_n_individuals[(i_seg-1)];
	d_mean_sum = 0;
	for(i_incr=0; i_incr<i_seg; i_incr++)
	  d_mean_sum += (i_n_individuals[i_incr] * d_standardised_x[i_incr]);
	d_mean_sum += (d_population_against_threshold * d_standardised_x[i_seg]);
	arg_d_revised_Z_upper_limit[i_seg] = (d_mean_sum / d_actual_population);
	d_var_sum = 0;
	for(i_incr=0; i_incr<i_seg; i_incr++)
	  d_var_sum += (i_n_individuals[i_incr] * ((d_standardised_x[i_incr] - arg_d_revised_Z_upper_limit[i_seg]) * (d_standardised_x[i_incr] - arg_d_revised_Z_upper_limit[i_seg])));
	d_var_sum += (d_population_against_threshold * ((d_standardised_x[i_seg] - arg_d_revised_Z_upper_limit[i_seg]) * (d_standardised_x[i_seg] - arg_d_revised_Z_upper_limit[i_seg])));
	arg_d_revised_Vp_upper_limit[i_seg] = (d_var_sum / d_actual_population);
	} // end for i_seg
  arg_d_revised_Z_upper_limit[(arg_i_n_segments-1)] = 0;
  arg_d_revised_Vp_upper_limit[(arg_i_n_segments-1)] = 1;

// TEMP - write the revised limits to file
/*
TempOutFileStream << "std_x n_individuals rev_Z_upperlimit rev_Vp_upperlimit" << "\n";
for(i_seg=0; i_seg<arg_i_n_segments; i_seg++)
  {
  TempOutFileStream << d_standardised_x[i_seg] << " ";
  TempOutFileStream << i_n_individuals[i_seg] << " ";
  TempOutFileStream << arg_d_revised_Z_upper_limit[i_seg] << " ";
  TempOutFileStream << arg_d_revised_Vp_upper_limit[i_seg] << "\n";
  } // end for i_seg
TempOutFileStream.close();
*/
// END TEMP

  // Free local arrays
  FreeDouble1DArray(d_standardised_x);
  FreeInt1DArray(i_n_individuals);

} // end Calculate_upper_threshold_correction

//------------------------------------------------------------------------------
// a simple function to round to the nearest integer
int round_double(double arg_d_number)
{
  int i_floor;
  double d_diff;

  i_floor = (int) (arg_d_number);
  d_diff = arg_d_number - i_floor;
  if(d_diff < 0.5)
	return(i_floor);
  else
	return(i_floor + 1);

} // end round_double

//------------------------------------------------------------------------------
// DISPERSAL MODE 1 CORRECTION - UNIVERSAL AVERAGING WITHIN CATEGORIES //
void Implement_dispersal_mode_1_correction(int arg_n_rows,
										   int arg_n_cols,
										   int arg_i_time,
										   int arg_i_n_dispersal_categories,
										   int** arg_i_dispersal_category,
										   short** arg_s_time_last_present,
										   double** arg_d_current_z,     // Gets changed here
										   double** arg_d_current_Vp)    // Gets changed here
{
  // Declare local variables
  int i_env, i_dis, i_row, i_col;
  int i_ncats_env = 100;
  double d_Z = 0;
  double d_Vp = 0;
  double d_x;
  double d_prob_so_far;
  double d_prop_above;
  double d_cat_prob;
  double d_variance_sum;
  double d_new_Vp;
  double d_data;
  double d_check_summed_prob;
  float f_Z;
  float f_Vp;

  // pointers to local arrays
  double** d_mean_dispersal;
  double** d_summed_prob;
  double* d_min_Z;
  double* d_max_Z;
  double* d_Vp_at_min_Z;
  double* d_Vp_at_max_Z;
  double* d_min_Vp;
  double* d_max_Vp;
  double* d_min_env;
  double* d_max_env;
  double* d_env_cat_increment;
  float f_data;

  // allocate local arrays
  d_mean_dispersal=AllocateDouble2DArray(arg_i_n_dispersal_categories,
										 3);
  d_summed_prob=AllocateDouble2DArray(i_ncats_env,
									  arg_i_n_dispersal_categories);
  d_min_Z=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_max_Z=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_Vp_at_min_Z=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_Vp_at_max_Z=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_min_Vp=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_max_Vp=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_min_env=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_max_env=AllocateDouble1DArray(arg_i_n_dispersal_categories);
  d_env_cat_increment=AllocateDouble1DArray(arg_i_n_dispersal_categories);

  // initialise the arrays
  for(i_dis=0; i_dis<arg_i_n_dispersal_categories; i_dis++)
	{
	d_mean_dispersal[i_dis][0]=0;
	d_mean_dispersal[i_dis][1]=0;
	d_mean_dispersal[i_dis][2]=0;
	d_min_Z[i_dis]=100000000;
	d_max_Z[i_dis]=-9999;
	d_min_Vp[i_dis]=100000000;
	d_max_Vp[i_dis]=-9999;
	for(i_env=0; i_env<i_ncats_env; i_env++)
	  d_summed_prob[i_env][i_dis]=0;
	} // end for i_dis

  // Read through the grid the first time, to calculate the average Z for
  // each category
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  if(arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1))
		{
		i_dis=arg_i_dispersal_category[i_row][i_col];
		i_dis -= 1;
		f_Z = (float)arg_d_current_z[i_row][i_col];
		f_Vp = (float)arg_d_current_Vp[i_row][i_col];
		d_Z = 0;
		d_Vp = 0;
		d_Z += (f_Z);
		d_Vp += (f_Vp);
		d_mean_dispersal[i_dis][0] += 1;
		d_mean_dispersal[i_dis][1] += d_Z;
		if(d_Z>d_max_Z[i_dis])
		  {
		  d_max_Z[i_dis] = d_Z;
		  d_Vp_at_max_Z[i_dis] = d_Vp;
		  } // end if arg_d_current_z[i_row][i_col]>d_max_Z[i_dis]
		if(d_Z<d_min_Z[i_dis])
		  {
		  d_min_Z[i_dis] = d_Z;
		  d_Vp_at_min_Z[i_dis] = d_Vp;
		  }
		if(d_Vp > d_max_Vp[i_dis])
		  d_max_Vp[i_dis] = d_Vp;
		if(d_Vp < d_min_Vp[i_dis])
		  d_min_Vp[i_dis] = d_Vp;
		} // end if s_time_last_present[i_row][i_col] == (i_time + 1)
	  } // end for i_col
	} // end for i_row
  // Now write these average Z for each dispersal category
  for(i_dis=0; i_dis<arg_i_n_dispersal_categories; i_dis++)
	{
	if(d_mean_dispersal[i_dis][0] > 0)
	  d_mean_dispersal[i_dis][1] = d_mean_dispersal[i_dis][1] / d_mean_dispersal[i_dis][0];
	d_min_env[i_dis] = d_min_Z[i_dis] - (d_Vp_at_min_Z[i_dis] * 8);
	d_max_env[i_dis] = d_max_Z[i_dis] + (d_Vp_at_max_Z[i_dis] * 8);
	d_env_cat_increment[i_dis] = (d_max_env[i_dis] - d_min_env[i_dis]) / (double) (i_ncats_env);
	d_data = d_min_env[i_dis]; //TEMP
	d_data = d_max_env[i_dis]; //TEMP
	d_data = d_env_cat_increment[i_dis]; //TEMP
	d_data += 1; //TEMP
	} // end for i_dis

  // Now loop through the dispersal categories. Use the Z & Vp for each active
  // grid cell in that dispersal category to calculate the new variance
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  if(arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1))
		{
		i_dis=arg_i_dispersal_category[i_row][i_col];
		i_dis -= 1;
		// Check whether we need to determine Vp for this dispersal category
		if(d_mean_dispersal[i_dis][2] <= 0)
		  {
		  f_Z = (float)arg_d_current_z[i_row][i_col];
		  f_Vp = (float)arg_d_current_Vp[i_row][i_col];
		  d_Z = 0;
		  d_Vp = 0;
		  d_Z += (f_Z);
		  d_Vp += (f_Vp);
		  d_x = d_min_env[i_dis];
		  d_prob_so_far = 0;
		  for(i_env=0; i_env<i_ncats_env; i_env++)
			{
			d_prop_above = calculate_above_limit_proportion(d_x,
															d_Z,
															d_Vp);
			d_cat_prob = ((1-d_prop_above)-d_prob_so_far);
			d_summed_prob[i_env][i_dis] += d_cat_prob;
			//f_data = (float)d_summed_prob[i_env][i_dis]; //TEMP
			d_prob_so_far += d_cat_prob;
			d_x += d_env_cat_increment[i_dis];
			} // end for i_env
		  //f_data += 1; //TEMP
		  } // end if d_mean_dispersal[i_dis][2] <= 0
		} // end if arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1)
	  } // end for i_col
	} // end for i_row

  // Now calculate the variance within each dispersal category
  for(i_dis=0; i_dis<arg_i_n_dispersal_categories; i_dis++)
	{
	if(d_mean_dispersal[i_dis][2] <= 0 && d_mean_dispersal[i_dis][0] > 0)
	  {
	  d_variance_sum = 0;
	  d_check_summed_prob = 0;
	  d_x = d_min_env[i_dis] - (d_env_cat_increment[i_dis]/2); // correct the x-value to sit in the middle of the increment
	  for(i_env=0; i_env<i_ncats_env; i_env++)
		{
		d_variance_sum += (((d_x - d_mean_dispersal[i_dis][1])*(d_x - d_mean_dispersal[i_dis][1])) * d_summed_prob[i_env][i_dis]);
		d_x += d_env_cat_increment[i_dis];
		d_check_summed_prob += d_summed_prob[i_env][i_dis];
		} // end for i_dis
	  d_new_Vp = d_variance_sum / d_check_summed_prob;
	  d_mean_dispersal[i_dis][2] = d_new_Vp;
	  } // end if d_mean_dispersal[i_dis][2] <= 0
	} // end for i_dis

  // Now put the right variance in each grid cell
  // Write the averages back to the Z & Vp arrays
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  if(arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1))
		{
		i_dis=arg_i_dispersal_category[i_row][i_col];
        i_dis -= 1;
		arg_d_current_z[i_row][i_col]  = d_mean_dispersal[i_dis][1];
		arg_d_current_Vp[i_row][i_col] = d_mean_dispersal[i_dis][2];
		} // end if arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1)
	  } // end for i_col
	} // end for i_row

  // Free the local array
  FreeDouble2DArray(d_mean_dispersal);
  FreeDouble2DArray(d_summed_prob);
  FreeDouble1DArray(d_min_Z);
  FreeDouble1DArray(d_max_Z);
  FreeDouble1DArray(d_Vp_at_min_Z);
  FreeDouble1DArray(d_Vp_at_max_Z);
  FreeDouble1DArray(d_min_Vp);
  FreeDouble1DArray(d_max_Vp);
  FreeDouble1DArray(d_min_env);
  FreeDouble1DArray(d_max_env);
  FreeDouble1DArray(d_env_cat_increment);

} // end Implement_dispersal_mode_1_correction
//------------------------------------------------------------------------------
// DISPERSAL MODE 2 CORRECTION - neighbourhood averaging by relative weight //
void Implement_dispersal_mode_2_correction(int arg_n_rows,
										   int arg_n_cols,
										   int arg_n_rel,
										   int arg_i_time,
										   short** arg_s_time_last_present,
										   int** arg_i_rel_row_col,
										   double* arg_d_rel_prob,
										   double** arg_d_current_z,     // Gets changed here
										   double** arg_d_current_Vp)    // Gets changed here
{
  // declare local variables
  int i_row, i_col, i_rel, i_cat;
  int i_ncats_env = 100;
  int i_relrow;
  int i_relcol;
  double d_sum_weight;
  double d_sum_Z;
  double d_sum_Vp;
  double d_max_Z, d_min_Z;
  double d_Vp_at_max_Z, d_Vp_at_min_Z;
  double d_min_env, d_max_env;
  double d_max_Vp;
  double d_min_Vp;
  double d_new_Z;
  double d_env_cat_increment;
  double d_Z;
  double d_Vp;
  double d_x;
  double d_prob_so_far;
  double d_prop_above;
  double d_cat_prob;
  double d_variance_sum;
  double d_new_Vp;
  double d_check_summed_prob;

  // allocate a local array to catch the Z & Vp of valid relative cells
  double** d_rel_Z_Vp;
  double* d_summed_prob;
  d_rel_Z_Vp = AllocateDouble2DArray(arg_n_rel,
									 2);
  d_summed_prob = AllocateDouble1DArray(i_ncats_env);

// Loop through all the cells on the grid, calculate and make the correction
for(i_row=0; i_row<arg_n_rows; i_row++)
  {
  for(i_col=0; i_col<arg_n_cols; i_col++)
	{
	if(arg_s_time_last_present[i_row][i_col] == (arg_i_time + 1))
	  {
	  d_sum_weight = 0;
	  d_sum_Z = 0;
	  d_sum_Vp = 0;
	  d_max_Z = -9999;
	  d_min_Z = 100000000;
	  d_max_Vp = -9999;
	  d_min_Vp = 100000000;
	  // initialise the array to catch Z & Vp
	  for(i_rel=0; i_rel<arg_n_rel; i_rel++)
		{
		d_rel_Z_Vp[i_rel][0] = -9999;
		d_rel_Z_Vp[i_rel][1] = -9999;
		} //end for i_rel
	  for(i_rel=0; i_rel<arg_n_rel; i_rel++)
		{
		i_relrow = i_row + arg_i_rel_row_col[i_rel][0];
		i_relcol = i_col + arg_i_rel_row_col[i_rel][1];
		// check if this is in the domain
		if(i_relrow > -1 && i_relrow < arg_n_rows)
		  {
		  if(i_relcol > -1 && i_relcol < arg_n_cols)
			{
			if(arg_s_time_last_present[i_relrow][i_relcol] == (arg_i_time + 1))
			  {
			  d_sum_weight += arg_d_rel_prob[i_rel];
			  d_sum_Z += (arg_d_current_z[i_relrow][i_relcol] * arg_d_rel_prob[i_rel]);
			  d_rel_Z_Vp[i_rel][0] = arg_d_current_z[i_relrow][i_relcol];
			  d_rel_Z_Vp[i_rel][1] = arg_d_current_Vp[i_relrow][i_relcol];
			  if(arg_d_current_z[i_relrow][i_relcol]>d_max_Z)
				{
				d_max_Z = arg_d_current_z[i_relrow][i_relcol];
				d_Vp_at_max_Z = arg_d_current_Vp[i_relrow][i_relcol];
				} // end if arg_d_current_z[i_relrow][i_relcol]>d_max_Z
			  if(arg_d_current_z[i_relrow][i_relcol]<d_min_Z)
				{
				d_min_Z = arg_d_current_z[i_relrow][i_relcol];
				d_Vp_at_min_Z = arg_d_current_Vp[i_relrow][i_relcol];
				}
			  if(arg_d_current_Vp[i_relrow][i_relcol] > d_max_Vp)
				d_max_Vp = arg_d_current_Vp[i_relrow][i_relcol];
			  if(arg_d_current_Vp[i_relrow][i_relcol] < d_min_Vp)
				d_min_Vp = arg_d_current_Vp[i_relrow][i_relcol];
			  } // end if s_time_last_present[i_relrow][i_relcol] == (i_time + 1)
			} // end if i_relcol > -1 && i_relcol < n_cols
		  } // end if i_relrow > -1 && i_relrow <n_rows
		} // end for i_rel
	  // calculate the mean Z and add it to the array, then work out the average
	  // variance (Vp) as well
	  if(d_sum_weight>0)
		{
		d_new_Z = d_sum_Z / d_sum_weight;
		arg_d_current_z[i_row][i_col] = d_new_Z;
		// Calculate the new Vp
		d_min_env = d_min_Z - (d_Vp_at_min_Z * 8);
		d_max_env = d_max_Z + (d_Vp_at_max_Z * 8);
		d_env_cat_increment = (d_max_env - d_min_env) / (double) (i_ncats_env);
		// initialise the probability summed catcher array
		for(i_cat=0; i_cat<i_ncats_env; i_cat++)
		  d_summed_prob[i_cat] = 0;
		// loop through the relative cells and add the relevent normal distributions
		// to the summed catcher array
		for(i_rel=0; i_rel<arg_n_rel; i_rel++)
		  {
		  if(d_rel_Z_Vp[i_rel][0] > -9999)
			{
			d_Z = d_rel_Z_Vp[i_rel][0];
			d_Vp = d_rel_Z_Vp[i_rel][1];
			d_x = d_min_env;
			d_prob_so_far = 0;
			for(i_cat=0; i_cat<i_ncats_env; i_cat++)
			  {
			  d_prop_above = calculate_above_limit_proportion(d_x,
															  d_Z,
															  d_Vp);
			  d_cat_prob = ((1-d_prop_above)-d_prob_so_far);
			  d_summed_prob[i_cat] += (d_cat_prob * arg_d_rel_prob[i_rel]);
			  d_prob_so_far += d_cat_prob;
			  d_x += d_env_cat_increment;
			  } // end for i_cat
			} // end if d_rel_Z_Vp[i_rel][0] > -9999
		  } // end for i_rel
		// Now use the summed weighted probabilities to calculate the new Vp
		d_variance_sum = 0;
		d_x = d_min_env - (d_env_cat_increment/2); // correct the x-value to sit in the middle of the increment
		d_check_summed_prob = 0;
		for(i_cat=0; i_cat<i_ncats_env; i_cat++)
		  {
		  d_variance_sum += (((d_x - d_new_Z)*(d_x - d_new_Z))*d_summed_prob[i_cat]);
		  d_check_summed_prob += d_summed_prob[i_cat];
		  d_x += d_env_cat_increment;
		  } // end for i_cat
		d_new_Vp = d_variance_sum / d_check_summed_prob;
		arg_d_current_Vp[i_row][i_col] = d_new_Vp;
		} // end if d_sum_weight>0
	  } // end if s_time_last_present[i_row][i_col] == (i_time + 1)
	} // end for i_col
  } // end for i_row

  // Free the allocated array
  FreeDouble2DArray(d_rel_Z_Vp);
  FreeDouble1DArray(d_summed_prob);

} // end Implement_dispersal_mode_2_correction
//------------------------------------------------------------------------------
// Assess persistence and adaptation in occupied grid cells
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
									   double* arg_d_nat_sel_slope)
{
  int i_row, i_col, i_env;
  int b_still_present;
  int i_array_index;
  int i_temp; //TEMP
  short int s_present_time;
  short int s_occurrence;
  short int s_zero = 0;
  short int s_new;
  double d_Env1;
  double d_S;
  double d_prop_above, d_prop_below;
  double d_relative_x;
  double d_relative_standard_x;
  double d_i;
  double d_delta_Z;
  double d_Vp_new;
  double d_z_diff;
  double d_i_fitness;

  // Pointers to GridCells
  TGridCell_ptr FromCell;

  // pointers to arrays on GridCells
  float* f_env_array;
  double* d_Z_env_array;
  double* d_Vp_env_array;

  // Loop through all the grid cells
  s_present_time = (short) (arg_i_time+1);
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  FromCell = Raft->GetCell(i_row,
							   i_col);
	  if(FromCell != NULL) // Check if this cell is land
		{
		// check if the species occurs in this cell
		s_occurrence = FromCell->GetOccurrence();
		if(s_occurrence>0)
		  {
		  b_still_present = 1;
		  f_env_array = FromCell->GetEnvData();
		  d_Z_env_array = FromCell->GetCurrentZ();
		  d_Vp_env_array = FromCell->GetCurrentVp();
		  for(i_env=0;i_env<arg_n_environments;i_env++)
			{
			//check if there is the potential for adaptation to this environment
			if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			  {
			  // Check that Vp (variance of attribite) is >0
			  if(d_Vp_env_array[i_env] > 0.000001)
				{
				// Determine which direction we're allowing adaptation, then
				// run the process
				if(arg_b_env_high_adaptation[i_env]>0) // adaptation at high values
				  {
			  /////// NEW - account for natural selection against higher values
				  if(arg_d_nat_sel_slope[i_env]>0)
				  {
				  if(d_Z_env_array[i_env]>arg_f_env_high_fund_limit[i_env])
					{
					d_z_diff = d_Z_env_array[i_env]-arg_f_env_high_fund_limit[i_env];
					d_i_fitness = d_z_diff * d_z_diff * arg_d_nat_sel_slope[i_env];
					d_delta_Z = d_i_fitness * arg_d_spp_H2[i_env] * sqrt(d_Vp_env_array[i_env]);
					if(d_delta_Z<d_z_diff)
					   d_Z_env_array[i_env] = d_Z_env_array[i_env] - d_delta_Z;
					else
					   d_Z_env_array[i_env] = arg_f_env_high_fund_limit[i_env];
					} // end if d_Z_env_array[i_env]>arg_f_env_high_fund_limit[i_env]
				  } // end if arg_d_nat_sel_slope[i_env]>0	
			 /////// END NEW
				  // only run adaptation if we're above the min tolerance limit & below the upper adaptability threshold
				  if(arg_f_env_low_fund_limit[i_env]<f_env_array[i_env] && arg_d_adapt_limit[i_env]>f_env_array[i_env])
					{
					// Calculate survival percentage (S)
					d_Env1=(double) (f_env_array[i_env]);
					d_S = calculate_survival_percentage(d_Env1,
														d_Z_env_array[i_env],
														d_Vp_env_array[i_env]);
					// Have a critical survival level in order to persist
					if(d_S>arg_d_spp_survival_min)
					  {
					  // But don't assess adaptation unless the survival percentage is
					  // notably less than 100%
					  if(d_S<(100-arg_d_spp_survival_min))
						{
						// Calculate percentage above the limit (d_spp_limit)
						d_prop_above = calculate_above_limit_proportion(arg_d_adapt_limit[i_env],
																		  d_Z_env_array[i_env],
																		  d_Vp_env_array[i_env]);
						// if the proportion of the population above the critical limit is
						// more than trivial (>0.01), calculate the actual d_Z and the
						// actual d_spp_Vp
						if(d_prop_above > 0.01)
						  {
						  // Extract actual mean (Z) and variance (Vp) for the upper-threshold
						  // distribution
						  d_relative_x = arg_d_adapt_limit[i_env] - d_Z_env_array[i_env];
						  d_relative_standard_x = d_relative_x / d_Vp_env_array[i_env];
						  i_array_index = (int) ((d_relative_standard_x * arg_d_increment_per_index) + arg_i_index_for_mean);
						  if(i_array_index > -1 && i_array_index < arg_i_n_segments)
							{
							//d_Z_env_array[i_env] = d_Z_env_array[i_env] + (d_Vp_env_array[i_env] * arg_d_revised_Z_upper_limit[i_array_index]); // REMOVE UPWARD DRIFT IN Z
							d_Vp_new = d_Vp_env_array[i_env] * arg_d_revised_Vp_upper_limit[i_array_index];
							} // end if i_array_index > -1 && i_array_index < i_n_segments
						  } // end if d_prop_above > 0.0001
						else
						  d_Vp_new = d_Vp_env_array[i_env];
						// Calculate i
						d_i = 2.2014 - (0.0488*d_S)+(0.00058*(d_S*d_S))-(0.0000029*(d_S*d_S*d_S));
						// Calculate delta Z
						d_delta_Z = d_i * arg_d_spp_H2[i_env] * sqrt(d_Vp_new);
						// Calculate the new Z & write it to the Z array						
						if(d_delta_Z > 0)
						  {
						  d_Z_env_array[i_env] = d_Z_env_array[i_env] + d_delta_Z;	 // NEW, MOVED WITHIN THE CHECK ON DELTA_Z BEING POSITIVE 
						  d_Vp_env_array[i_env] = d_Vp_new;
//TEMP
if(d_Vp_new>1.106)
  i_temp=1;
//TEMP
						  }	  
						} // end if d_S<(100-arg_d_spp_survival_min)
					  } // end if d_S>d_spp_survival_min
					else
					  b_still_present = b_still_present * 0;
					} // end if f_env_high_fund_limit[i_env]>f_env_array[i_env]
				  else
					b_still_present = b_still_present * 0;
				  } // end if arg_b_env_high_adaptation>0
				else // adaptation at low values
				  {
			 /////// NEW - account for natural selection against lower values
				if(arg_d_nat_sel_slope[i_env]>0)
				  {   
				  if(d_Z_env_array[i_env]<arg_f_env_low_fund_limit[i_env])
					{
					d_z_diff = arg_f_env_low_fund_limit[i_env]-d_Z_env_array[i_env];
					d_i_fitness = d_z_diff * d_z_diff * arg_d_nat_sel_slope[i_env];
					d_delta_Z = d_i_fitness * arg_d_spp_H2[i_env] * sqrt(d_Vp_env_array[i_env]);
					if(d_delta_Z<d_z_diff)
					   d_Z_env_array[i_env] = d_Z_env_array[i_env] + d_delta_Z;
					else
					   d_Z_env_array[i_env] = arg_f_env_low_fund_limit[i_env];
					} // end if d_Z_env_array[i_env]>arg_f_env_high_fund_limit[i_env]
				  } // end if arg_d_nat_sel_slope[i_env]>0					
			 /////// END NEW
				  // only run adaptation if we're below the max fundamental limit & above the low adaptability limit
				  if(arg_f_env_high_fund_limit[i_env]>f_env_array[i_env] && arg_d_adapt_limit[i_env]<f_env_array[i_env])
					{
					// Calculate survival percentage (S)
					d_Env1=(double) (f_env_array[i_env]);
					d_S = calculate_survival_percentage_LowLimit(d_Env1,
																 d_Z_env_array[i_env],
																 d_Vp_env_array[i_env]);
					// Have a critical survival level in order to persist
					if(d_S>arg_d_spp_survival_min)
					  {
                      // But don't assess adaptation unless the survival percentage is
					  // notably less than 100%
					  if(d_S<(100-arg_d_spp_survival_min))
						{
						// Calculate percentage above the limit (d_spp_limit)
						d_prop_below = calculate_below_limit_proportion(arg_d_adapt_limit[i_env],
																		  d_Z_env_array[i_env],
																		  d_Vp_env_array[i_env]);
						// if the proportion of the population above the critical limit is
						// more than trivial (>0.01), calculate the actual d_Z and the
						// actual d_spp_Vp
						if(d_prop_below > 0.01)
						  {
						  // Extract actual mean (Z) and variance (Vp) for the upper-threshold
						  // distribution
						  d_relative_x = d_Z_env_array[i_env] - arg_d_adapt_limit[i_env];
						  d_relative_standard_x = d_relative_x / d_Vp_env_array[i_env];
						  i_array_index = (int) ((d_relative_standard_x * arg_d_increment_per_index) + arg_i_index_for_mean);
						  if(i_array_index > -1 && i_array_index < arg_i_n_segments)
							{
							//d_Z_env_array[i_env] = d_Z_env_array[i_env] - (d_Vp_env_array[i_env] * arg_d_revised_Z_upper_limit[i_array_index]);
							d_Vp_new = d_Vp_env_array[i_env] * arg_d_revised_Vp_upper_limit[i_array_index];
							} // end if i_array_index > -1 && i_array_index < i_n_segments
						  } // end if d_prop_below > 0.0001
						else
						  d_Vp_new = d_Vp_env_array[i_env];
						// Calculate i
						d_i = 2.2014 - (0.0488*d_S)+(0.00058*(d_S*d_S))-(0.0000029*(d_S*d_S*d_S));
						// Calculate delta Z
						d_delta_Z = d_i * arg_d_spp_H2[i_env] * sqrt(d_Vp_new);
						// Calculate the new Z & write it to the Z array						
						if(d_delta_Z > 0)
						  {
						  d_Z_env_array[i_env] = d_Z_env_array[i_env] - d_delta_Z; // NEWLY MOVED HERE TO AVOID CHANGE IN THE THRESHOLD THE WRONG WAY
						  d_Vp_env_array[i_env] = d_Vp_new;
						  } // end if d_delta_Z > 0
						} // end if d_S<(100-arg_d_spp_survival_min)
					  } // end if d_S>d_spp_survival_min
					else
					  b_still_present = b_still_present * 0;
					} // end if f_env_low_fund_limit[i_env]<f_env_array[i_env]
				  else
					 b_still_present = b_still_present * 0;
				  } // end else b_env_low_adaptation[i_env]>0
				} // end if d_Vp_env_array[i_env] > 0
			  else
				b_still_present = b_still_present * 0;
			  } // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
			else // no adaptation
			  {
			  if(arg_f_env_low_fund_limit[i_env]<f_env_array[i_env] && arg_f_env_high_fund_limit[i_env]>f_env_array[i_env])
				b_still_present = b_still_present * 1;
			  else
				b_still_present = b_still_present * 0;
			  } // end else b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
			} // end for i_env
		  if(b_still_present>0) // if the species has persisted in this grid cell
			{
			// increment the time last present record
			FromCell->SetTimeLastPresent(s_present_time);
			} // end b_still_present>0
		  else
			FromCell->SetOccurrence(s_zero);
		  } // end if s_occurrence>0
		} // end if FromCell != NULL
	  } // end for i_col
	} // end for i_row

} // end Assess_Persistence_and_Adaptation

//------------------------------------------------------------------------------
// Implement dispersal to unnoccupied grid cells
// Dispersal success is based on the environmental conditions in each potential
// destination cell being suitable for organisms arriving from a given source cell.
// Dispersal probabilities are as specified in the relative neighbourhood file.
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
											 double* arg_d_adapt_limit)
{
  int i_row, i_col, i_env, i_rel;
  int i_ncats_env = 100;
  int b_environment_suitable;
  int i_relrow, i_relcol;
  short int s_occurrence;
  short int s_one = 1;
  float f_random;
  float f_min = 0;
  float f_max = 1;
  double d_prob_sum;
  double d_aggregate_dispersal_prob;
  double d_dispersal_complement;
  double d_aggregate_dispersal_complement;

  // Pointers to GridCells
  TGridCell_ptr ToCell;
  TGridCell_ptr FromCell;

  // pointers to arrays on GridCells
  float* f_env_array;
  double* d_Z_env_array;
  double* d_Vp_env_array;
  double* d_From_Z_env_array;
  double* d_From_Vp_env_array;

// allocate a local array to catch the Z & Vp of valid relative cells
  double** d_rel_Z;
  double** d_rel_Vp;
  double* d_summed_prob;
  d_rel_Z = AllocateDouble2DArray(arg_n_rel,
								  arg_n_environments);
  d_rel_Vp = AllocateDouble2DArray(arg_n_rel,
								  arg_n_environments);
  d_summed_prob = AllocateDouble1DArray(i_ncats_env);

  // Loop through the grid to assess dispersal
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  ToCell = Raft->GetCell(i_row,
							 i_col);
	  if(ToCell != NULL) // Check if this cell is land
		{
		// check if the species occurs in this cell
		s_occurrence = ToCell->GetOccurrence();
		//if(s_occurrence<1) // the species is absent  //OLD
		if(s_occurrence == 0) // the species is absent //NEW 17 July 2016
		  {
		  b_environment_suitable = 1;
		  f_env_array = ToCell->GetEnvData();
		  d_Z_env_array = ToCell->GetCurrentZ();
		  d_Vp_env_array = ToCell->GetCurrentVp();
		  for(i_env=0;i_env<arg_n_environments;i_env++)
			{
			//check if there is the potential for adaptation to this environment
			if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			  {
			  if(arg_b_env_low_adaptation[i_env]>0) // adaptation at low values
				{
				if(arg_f_env_high_fund_limit[i_env]>f_env_array[i_env] && arg_d_adapt_limit[i_env]<f_env_array[i_env])
				  b_environment_suitable = b_environment_suitable * 1;
				else
				  b_environment_suitable = b_environment_suitable * 0;
				} // end if b_env_low_adaptation[i_env]>0
			  else // adaptation at high values
				{
				if(arg_f_env_low_fund_limit[i_env]<f_env_array[i_env] && arg_d_adapt_limit[i_env]>f_env_array[i_env])
				  b_environment_suitable = b_environment_suitable * 1;
				else
				  b_environment_suitable = b_environment_suitable * 0;
				} // end else b_env_low_adaptation[i_env]>0
			  } // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
			else // no adaptation on this environment
			  {
			  if(arg_f_env_low_fund_limit[i_env]<f_env_array[i_env] && arg_f_env_high_fund_limit[i_env]>f_env_array[i_env])
				b_environment_suitable = b_environment_suitable * 1;
			  else
				b_environment_suitable = b_environment_suitable * 0;
			  } // end else b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
			} // end for i_env
		  // Now, if the environment in this grid cell is suitable, look in
		  // the dispersal neighbourhood & determine the probability of the
		  // species colonising this cell.
		  if(b_environment_suitable>0)
			{
			// Initialise the catching variables & arrays
			d_prob_sum = 0;
			// Now loop through the relative neighbourhood around this cell,
			// catch & aggregate the dispersal probability & the adaptive
			// variables
			for(i_rel=0; i_rel<arg_n_rel; i_rel++)
			  {
			  i_relrow = i_row + arg_i_rel_row_col[i_rel][0];
			  i_relcol = i_col + arg_i_rel_row_col[i_rel][1];
			  // ensure the relative cell is still on the raft
			  if(i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols)
				{
				FromCell = Raft->GetCell(i_relrow,
										 i_relcol);
				if(FromCell != NULL) // Check if this cell is land
				  {
				  s_occurrence = FromCell->GetOccurrence();
				  if(s_occurrence>0) // the species is present
					{
					d_dispersal_complement = 1 - arg_d_rel_prob[i_rel];  //NEW 9 Apr 15
					if(d_prob_sum>0) // we already have                   //NEW 9 Apr 15
					  d_aggregate_dispersal_complement = d_aggregate_dispersal_complement * d_dispersal_complement; //NEW 9 Apr 15
					else                                                 //NEW 9 Apr 15
					  d_aggregate_dispersal_complement = d_dispersal_complement; //NEW 9 Apr 15
					d_prob_sum += arg_d_rel_prob[i_rel];
					d_From_Z_env_array = FromCell->GetCurrentZ();
					d_From_Vp_env_array = FromCell->GetCurrentVp();
					for(i_env=0;i_env<arg_n_environments;i_env++)
					  {
					  d_rel_Z[i_rel][i_env] = d_From_Z_env_array[i_env];
					  d_rel_Vp[i_rel][i_env] = d_From_Vp_env_array[i_env];
					  } // end for i_env
					} // end if s_occurrence>0
				  else
					{
					for(i_env=0;i_env<arg_n_environments;i_env++)
					  {
					  d_rel_Z[i_rel][i_env] = -9999;
					  d_rel_Vp[i_rel][i_env] = -9999;
					  } // end for i_env
					}// end else s_occurrence>0
				  } // end if FromCell != NULL
				else
				  {
				  for(i_env=0;i_env<arg_n_environments;i_env++)
					{
					d_rel_Z[i_rel][i_env] = -9999;
					d_rel_Vp[i_rel][i_env] = -9999;
					} // end for i_env
				  }// end else FromCell != NULL
				} // end if i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols
			  else //NEW 17 July 2016 - for when part of the relative neighbourhood is outside the grid
				{  												//NEW 17 July 2016
				for(i_env=0;i_env<arg_n_environments;i_env++)   //NEW 17 July 2016
				  {                                             //NEW 17 July 2016
				  d_rel_Z[i_rel][i_env] = -9999;   				//NEW 17 July 2016
				  d_rel_Vp[i_rel][i_env] = -9999;  				//NEW 17 July 2016
				  } // end for i_env                            //NEW 17 July 2016
				} //end else                                    //NEW 17 July 2016
			  }// end for i_rel
			// Now do a binomial trial to see if dispersal is successful, if
			// so, set this grid cell as occupied, determine the adaptive
			// parameters and set them to this grid cell.
			if(d_prob_sum>0)
			  {
			  d_aggregate_dispersal_prob = 1 - d_aggregate_dispersal_complement;  //NEW 9 Apr 15
			  // random throw
			  f_random = uniformRandomDeviate(f_min,
											  f_max);
			  if(d_aggregate_dispersal_prob>f_random)
				{
				// Mark the cell as newly colonised
				ToCell->SetNewlyColonised(s_one);
				// calculate mean and variance from contributing cells in the
				// neighbourhood & set these on this GridCell on the Raft
				Calculate_Newly_Colonised_Z_Vp_Mixture(Raft,
													   i_row,
													   i_col,
													   arg_n_rel,
													   arg_n_environments,
													   arg_b_env_low_adaptation,
													   arg_b_env_high_adaptation,
													   arg_f_env_low_fund_limit,
													   arg_f_env_high_fund_limit,
													   arg_d_rel_prob,
													   arg_d_adapt_limit,
													   d_rel_Z,
													   d_rel_Vp);

				} // end if d_prob_sum>f_random
			  } // end if d_prob_sum>0
			} // end if b_environment_suitable>0
		  } // end if s_occurrence<1
		} // end if ToCell != NULL
	  } // end for i_col
	} // end for i_row

} // end Implement_Dispersal_To_Unoccupied_Cells

//------------------------------------------------------------------------------
// A function to calculate the Z & Vp (attribute mean and it's variance) for
// the newly colonised grid cell, based on those grid cells contributing
// propagules.
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
									double** arg_d_rel_Vp)
{
  int i_row, i_col, i_env, i_rel, i_cat;
  double d_sum_Z;
  double d_sum_Vp;
  double d_max_Z;
  double d_min_Z;
  double d_max_Vp;
  double d_min_Vp;
  double d_new_Z, d_new_Vp;
  double d_min_env, d_max_env;
  double d_Vp_at_min_Z, d_Vp_at_max_Z;
  double d_env_cat_increment;
  double d_Z, d_Vp;
  double d_x;
  double d_prob_so_far;
  double d_prop_above;
  double d_cat_prob;
  double d_variance_sum;
  double d_check_summed_prob;

  // Pointers to GridCells
  TGridCell_ptr ToCell;

  //pointers to local arrays
  double* d_Z_env_array;
  double* d_Vp_env_array;

  // Get the focal cell
  ToCell = Raft->GetCell(arg_i_row,
						 arg_i_col);
  if(ToCell != NULL) // Double check this cell is land
	{
	d_Z_env_array = ToCell->GetCurrentZ();
	d_Vp_env_array = ToCell->GetCurrentVp();
	// Loop through each env for the focal cell. If it's an adaptive env, calculate
	// the mean (z) & variance (Vp)
	for(i_env=0;i_env<arg_n_environments;i_env++)
	  {
	  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
		{
		// initialise the catching data
		d_sum_Z = 0;
		d_sum_Vp = 0;
		d_max_Z = -9999;
		d_min_Z = 100000000;
		d_max_Vp = -9999;
		d_min_Vp = 100000000;
		for(i_rel=0; i_rel<arg_n_rel; i_rel++)
		  {
		  if(arg_d_rel_Z[i_rel][i_env] > -9999)
			{
			d_sum_Z += (arg_d_rel_Z[i_rel][i_env] * arg_d_rel_prob[i_rel]);
			if(arg_d_rel_Z[i_rel][i_env]>d_max_Z)
			  {
			  d_max_Z = arg_d_rel_Z[i_rel][i_env];
			  d_Vp_at_max_Z = arg_d_rel_Vp[i_rel][i_env];
			  } // end if arg_d_rel_Z[i_rel][i_env]>d_max_Z
			if(arg_d_rel_Z[i_rel][i_env]<d_min_Z)
			  {
			  d_min_Z = arg_d_rel_Z[i_rel][i_env];
			  d_Vp_at_min_Z = arg_d_rel_Vp[i_rel][i_env];
			  } // end if arg_d_rel_Z[i_rel][i_env]<d_min_Z
			if(arg_d_rel_Vp[i_rel][i_env] > d_max_Vp)
			  d_max_Vp = arg_d_rel_Vp[i_rel][i_env];
			if(arg_d_rel_Vp[i_rel][i_env] < d_min_Vp)
			  d_min_Vp = arg_d_rel_Vp[i_rel][i_env];
			} // end if arg_d_rel_Z[i_rel][i_env] > -9999
		  } // end for i_rel
		if((d_max_Z - d_min_Z < 0.0001) && (d_max_Vp - d_min_Vp < 0.0001))  // NEW 9 Apr 15
		  {                                                                     // NEW 9 Apr 15
		  d_Z_env_array[i_env] = d_max_Z; // NEW 9 Apr 15
		  d_Vp_env_array[i_env] = d_max_Vp; // NEW 9 Apr 15
		  } // end if d_max_Z - d_min_z < 0.0001) && (d_max_Vp - d_max_Vp < 0.0001)// NEW 9 Apr 15
		else    // NEW 9 Apr 15
		  {     // NEW 9 Apr 15
		  if(arg_d_prob_sum>0)
			{
			d_new_Z = d_sum_Z / arg_d_prob_sum;
			// Calculate the new Vp
			d_min_env = d_min_Z - (d_max_Vp * 4.5); //d_min_env = d_min_Z - (d_Vp_at_min_Z * 8); // <- changed from this
			// correct this to ensure variance doesn't drift beyond the lower limit 		   //NEW
			//if(arg_b_env_low_adaptation[i_env]>0 && d_min_env<arg_d_adapt_limit[i_env])      //NEW
			//  d_min_env = arg_d_adapt_limit[i_env];                                          //NEW
			d_max_env = d_max_Z + (d_max_Vp * 4.5);//d_max_env = d_max_Z + (d_Vp_at_max_Z * 8); // <- changed from this
			// correct this to ensure variance doesn't drift beyond the lower limit 		   //NEW
			//if(arg_b_env_high_adaptation[i_env]>0 && d_max_env>arg_d_adapt_limit[i_env])      //NEW
			//  d_max_env = arg_d_adapt_limit[i_env];                                          //NEW
			d_env_cat_increment = (d_max_env - d_min_env) / (double) (arg_i_ncats_env);
			// initialise the probability summed catcher array
			for(i_cat=0; i_cat<arg_i_ncats_env; i_cat++)
			  arg_d_summed_prob[i_cat] = 0;
			// loop through the relative cells and add the relevent normal distributions
			// to the summed catcher array
			for(i_rel=0; i_rel<arg_n_rel; i_rel++)
			  {
			  if(arg_d_rel_Vp[i_rel][i_env] > 0)
				{
				d_Z = arg_d_rel_Z[i_rel][i_env];
				d_Vp = arg_d_rel_Vp[i_rel][i_env];
				d_x = d_min_env;
				d_prob_so_far = 0;
				for(i_cat=0; i_cat<arg_i_ncats_env; i_cat++)
				  {
				  d_prop_above = calculate_above_limit_proportion(d_x,
																	d_Z,
																	d_Vp);
				  d_cat_prob = ((1-d_prop_above)-d_prob_so_far);
				  arg_d_summed_prob[i_cat] += (d_cat_prob * arg_d_rel_prob[i_rel]);
				  d_prob_so_far += d_cat_prob; /////////////////!!!!
				  d_x += d_env_cat_increment;
				  } // end for i_cat
				// put any left-over probability in the highest category          // NEW
				if(d_prob_so_far<1)                                               // NEW
				  arg_d_summed_prob[(arg_i_ncats_env-1)] += ((1-d_prob_so_far)*arg_d_rel_prob[i_rel]); // NEW
				} // end if d_rel_Z_Vp[i_rel][0] > -9999
			  } // end for i_rel
			// Now use the summed weighted probabilities to calculate the new Vp
			d_variance_sum = 0;
			d_x = d_min_env - (d_env_cat_increment/2); // correct the x-value to sit in the middle of the increment
			d_check_summed_prob = 0;
			for(i_cat=0; i_cat<arg_i_ncats_env; i_cat++)
			  {
			  d_variance_sum += (((d_x - d_new_Z)*(d_x - d_new_Z))*arg_d_summed_prob[i_cat]);
			  d_check_summed_prob += arg_d_summed_prob[i_cat];
			  d_x += d_env_cat_increment;
			  } // end for i_cat
			d_new_Vp = d_variance_sum / d_check_summed_prob;
			// ensure the new Z doesn't exceed the evolutionary limit
			if(arg_b_env_low_adaptation[i_env]>0)
			  {
			  if(d_new_Z<arg_d_adapt_limit[i_env])
				d_new_Z = arg_d_adapt_limit[i_env];
			  } // end if arg_b_env_low_adaptation[i_env]>0
			else
			  {
			  if(d_new_Z>arg_d_adapt_limit[i_env])
				d_new_Z = arg_d_adapt_limit[i_env];
			  } // end else arg_b_env_low_adaptation[i_env]>0
			// Put the calculated mean and variance in the right spot

			
			d_Z_env_array[i_env] = d_new_Z;
			d_Vp_env_array[i_env] = d_new_Vp;
			} // end if arg_d_prob_sum>0
		  } // end else      // NEW 9 Apr 15
		} // end if arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0
	  } // end for i_env
	} // end if ToCell != NULL

} // end Calculate_Newly_Colonised_Z_Vp

//------------------------------------------------------------------------------
// Apply adaptive parameter averaging
// Accounts for dispersal between occupied grid cells resulting in averaging of
// mean attributes and variance of attribute
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
											double* arg_d_adapt_limit)
{
  // declare local variables
  int i_row, i_col, i_rel, i_cat, i_env;
  int i_ncats_env = 100;
  int i_relrow;
  int i_relcol;
  short int s_occurrence;
  short int s_rel_occurrence;
  double d_sum_weight;
  double d_sum_Z;
  double d_sum_Vp;
  double d_max_Z, d_min_Z;
  double d_Vp_at_max_Z, d_Vp_at_min_Z;
  double d_min_env, d_max_env;
  double d_max_Vp;
  double d_min_Vp;
  double d_new_Z;
  double d_env_cat_increment;
  double d_Z;
  double d_Vp;
  double d_x;
  double d_prob_so_far;
  double d_prop_above;
  double d_cat_prob;
  double d_variance_sum;
  double d_new_Vp;
  double d_check_summed_prob;

  TGridCell_ptr ToCell;
  TGridCell_ptr FromCell;

  // pointers to arrays
  double* d_Z_env_array;
  double* d_Vp_env_array;
  double* d_next_Z_env_array;
  double* d_next_Vp_env_array;
  double* d_rel_Z_env_array;
  double* d_rel_Vp_env_array;

  // allocate a local array to catch the Z & Vp of valid relative cells
  double** d_rel_Z_Vp;
  double* d_summed_prob;
  d_rel_Z_Vp = AllocateDouble2DArray(arg_n_rel,
									 2);
  d_summed_prob = AllocateDouble1DArray(i_ncats_env);

// Loop through all the cells on the grid, calculate and make the correction
for(i_row=0; i_row<arg_n_rows; i_row++)
  {
  for(i_col=0; i_col<arg_n_cols; i_col++)
	{
	ToCell = Raft->GetCell(i_row,
						   i_col);
	if(ToCell != NULL) // Check if this cell is land
	  {
	  // check if the species occurs in this cell
	  s_occurrence = ToCell->GetOccurrence();
	  if(s_occurrence>0) // the species is present
		{
		// obtain the arrays for this grid cell
		d_Z_env_array = ToCell->GetCurrentZ();
		d_Vp_env_array = ToCell->GetCurrentVp();
		d_next_Z_env_array = ToCell->GetNextZ();
		d_next_Vp_env_array = ToCell->GetNextVp();
		for(i_env=0;i_env<arg_n_environments;i_env++)
		  {
		  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			{
			d_sum_weight = 0;
			d_sum_Z = 0;
			d_sum_Vp = 0;
			d_max_Z = -9999;
			d_min_Z = 100000000;
			d_max_Vp = -9999;
			d_min_Vp = 100000000;
			// initialise the array to catch Z & Vp
			for(i_rel=0; i_rel<arg_n_rel; i_rel++)
			  {
			  d_rel_Z_Vp[i_rel][0] = -9999;
			  d_rel_Z_Vp[i_rel][1] = -9999;
			  } //end for i_rel
			for(i_rel=0; i_rel<arg_n_rel; i_rel++)
			  {
			  i_relrow = i_row + arg_i_rel_row_col[i_rel][0];
			  i_relcol = i_col + arg_i_rel_row_col[i_rel][1];
			  // ensure the relative cell is still on the raft
			  if(i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols)
				{
				// check if this is in the domain
				FromCell = Raft->GetCell(i_relrow,
										 i_relcol);
				if(FromCell != NULL) // Check if this cell is land
				  {
				  // check if the species occurs in this cell
				  s_rel_occurrence = FromCell->GetOccurrence();
				  if(s_rel_occurrence>0) // the species is present
					{
					d_rel_Z_env_array = FromCell->GetCurrentZ();
					d_rel_Vp_env_array = FromCell->GetCurrentVp();
					d_sum_weight += arg_d_rel_prob[i_rel];
					d_sum_Z += (d_rel_Z_env_array[i_env] * arg_d_rel_prob[i_rel]);
					d_rel_Z_Vp[i_rel][0] = d_rel_Z_env_array[i_env];
					d_rel_Z_Vp[i_rel][1] = d_rel_Vp_env_array[i_env];
					if(d_rel_Z_env_array[i_env]>d_max_Z)
					  {
					  d_max_Z = d_rel_Z_env_array[i_env];
					  d_Vp_at_max_Z = d_rel_Vp_env_array[i_env];
					  } // end if d_Z_env_array[i_env]>d_max_Z
					if(d_rel_Z_env_array[i_env]<d_min_Z)
					  {
					  d_min_Z = d_rel_Z_env_array[i_env];
					  d_Vp_at_min_Z = d_rel_Vp_env_array[i_env];
					  }
					if(d_rel_Vp_env_array[i_env] > d_max_Vp)
					  d_max_Vp = d_rel_Vp_env_array[i_env];
					if(d_rel_Vp_env_array[i_env] < d_min_Vp)
					  d_min_Vp = d_rel_Vp_env_array[i_env];
					} // end if s_rel_occurrence>0
				  } // end if FromCell != NULL
				} // end if i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols
			  } // end for i_rel
			// Correct for the greater weighting the focal cell has in having a
			// larger resident population as opposed to a few successful dispersing
			// propagules
			d_sum_weight += arg_d_resident_weighting;
			d_sum_Z += (d_Z_env_array[i_env] * arg_d_resident_weighting);
			// calculate the mean Z and add it to the array, then work out the average
			// variance (Vp) as well
			if(d_sum_weight>0)
			  {
			  // A NEW CHECK TO SKIP THE VARIANCE AVERAGING WHERE APPROPRIATE - 8 Apr 15
			  if(((d_max_Z - d_min_Z) < 0.0001) && ((d_max_Vp - d_min_Vp) < 0.0001))
				{
				d_next_Z_env_array[i_env] = d_Z_env_array[i_env];
				d_next_Vp_env_array[i_env] = d_Vp_env_array[i_env];
				}
			  else
				{
			  d_new_Z = d_sum_Z / d_sum_weight;
			  // Calculate the new Vp
			  d_min_env = d_min_Z - (d_max_Vp * 4.5); //d_min_env = d_min_Z - (d_Vp_at_min_Z * 8); // <- changed from this
			  // correct this to ensure variance doesn't drift beyond the lower limit 		   //NEW
			  //if(arg_b_env_low_adaptation[i_env]>0 && d_min_env<arg_d_adapt_limit[i_env])      //NEW
			  //	d_min_env = arg_d_adapt_limit[i_env];                                          //NEW
			  d_max_env = d_max_Z + (d_max_Vp * 4.5);//d_max_env = d_max_Z + (d_Vp_at_max_Z * 8); // <- changed from this
			  // correct this to ensure variance doesn't drift beyond the lower limit 		   //NEW
			  // if(arg_b_env_high_adaptation[i_env]>0 && d_max_env>arg_d_adapt_limit[i_env])      //NEW
			  //	d_max_env = arg_d_adapt_limit[i_env];                                          //NEW
			  d_env_cat_increment = (d_max_env - d_min_env) / (double) (i_ncats_env);
			  // initialise the probability summed catcher array
			  for(i_cat=0; i_cat<i_ncats_env; i_cat++)
				d_summed_prob[i_cat] = 0;
			  // loop through the relative cells and add the relevent normal distributions
			  // to the summed catcher array
			  for(i_rel=0; i_rel<arg_n_rel; i_rel++)
				{
				if(d_rel_Z_Vp[i_rel][1] > 0) // only do this if the variance is >0
				  {
				  d_Z = d_rel_Z_Vp[i_rel][0];
				  d_Vp = d_rel_Z_Vp[i_rel][1];
				  d_x = d_min_env;
				  d_prob_so_far = 0;
				  for(i_cat=0; i_cat<i_ncats_env; i_cat++)
					{
					d_prop_above = calculate_above_limit_proportion(d_x,
																	d_Z,
																	d_Vp);
					d_cat_prob = ((1-d_prop_above)-d_prob_so_far);
					d_summed_prob[i_cat] += (d_cat_prob * arg_d_rel_prob[i_rel]);
					d_prob_so_far += d_cat_prob;
					d_x += d_env_cat_increment;
					} // end for i_cat
				  // put any left-over probability in the highest category          // NEW
				  if(d_prob_so_far<1)                                               // NEW
					d_summed_prob[(i_ncats_env-1)] += ((1-d_prob_so_far)*arg_d_rel_prob[i_rel]); // NEW
				  } // end if d_rel_Z_Vp[i_rel][0] > -9999
				} // end for i_rel
			  // Again, correct for the greater weighting the focal cell has in
			  // having a larger resident population as opposed to a few successful
			  // dispersing propagules
			  d_Z = d_Z_env_array[i_env];
			  d_Vp = d_Vp_env_array[i_env];
			  d_x = d_min_env;
			  d_prob_so_far = 0;
			  for(i_cat=0; i_cat<i_ncats_env; i_cat++)
				{
				d_prop_above = calculate_above_limit_proportion(d_x,
																d_Z,
																d_Vp);
				d_cat_prob = ((1-d_prop_above)-d_prob_so_far);
				d_summed_prob[i_cat] += (d_cat_prob * arg_d_resident_weighting);
				d_prob_so_far += d_cat_prob;
				d_x += d_env_cat_increment;
				} // end for i_cat
			  // put any left-over probability in the highest category          // NEW
			  if(d_prob_so_far<1)                                               // NEW
				d_summed_prob[(i_ncats_env-1)] += ((1-d_prob_so_far)*arg_d_resident_weighting); // NEW
			  // Now use the summed weighted probabilities to calculate the new Vp
			  d_variance_sum = 0;
			  d_x = d_min_env - (d_env_cat_increment/2); // correct the x-value to sit in the middle of the increment
			  d_check_summed_prob = 0;
			  for(i_cat=0; i_cat<i_ncats_env; i_cat++)
				{
				//TEMP REMOVE   d_variance_sum += (((d_x - d_new_Z)*(d_x - d_new_Z))*d_summed_prob[i_cat]);
				d_variance_sum += (((d_x - d_Z_env_array[i_env])*(d_x - d_Z_env_array[i_env]))*d_summed_prob[i_cat]);
				d_check_summed_prob += d_summed_prob[i_cat];
				d_x += d_env_cat_increment;
				} // end for i_cat
			  d_new_Vp = d_variance_sum / d_check_summed_prob;
			  // ensure the new Z doesn't exceed the evolutionary limit
			  if(arg_b_env_low_adaptation[i_env]>0)
				{
				if(d_new_Z<arg_d_adapt_limit[i_env])
				  d_new_Z = arg_d_adapt_limit[i_env];
				} // end if arg_b_env_low_adaptation[i_env]>0
			  else
				{
				if(d_new_Z>arg_d_adapt_limit[i_env])
				  d_new_Z = arg_d_adapt_limit[i_env];
				} // end else arg_b_env_low_adaptation[i_env]>0
			  // and finally record the new Z & Vp value
			  d_next_Z_env_array[i_env] = d_new_Z;
			  d_next_Vp_env_array[i_env] = d_new_Vp;
				} // end else d_max_Z == d_min_Z && d_max_Vp == d_min_Vp
			  } // end if d_sum_weight>0
			else
			  {
			  d_next_Z_env_array[i_env] = d_Z_env_array[i_env];
			  d_next_Vp_env_array[i_env] = d_Vp_env_array[i_env];
			  } // end else d_sum_weight>0
			} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
		  } // end for i_env
		} // end if s_occurrence>0
	  } // end if ToCell != NULL
	} // end for i_col
  } // end for i_row

// Now zip through the grid again and put the new Z & Vp values in the right spot
for(i_row=0; i_row<arg_n_rows; i_row++)
  {
  for(i_col=0; i_col<arg_n_cols; i_col++)
	{
	ToCell = Raft->GetCell(i_row,
						   i_col);
	if(ToCell != NULL) // Check if this cell is land
	  {
	  // check if the species occurs in this cell
	  s_occurrence = ToCell->GetOccurrence();
	  if(s_occurrence>0) // the species is present
		{
		// obtain the arrays for this grid cell
		d_Z_env_array = ToCell->GetCurrentZ();
		d_Vp_env_array = ToCell->GetCurrentVp();
		d_next_Z_env_array = ToCell->GetNextZ();
		d_next_Vp_env_array = ToCell->GetNextVp();
		for(i_env=0;i_env<arg_n_environments;i_env++)
		  {
		  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			{
			d_Z_env_array[i_env] = d_next_Z_env_array[i_env];
			d_Vp_env_array[i_env] = d_next_Vp_env_array[i_env];
			} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
		  } // end for i_env
		} // end if s_occurrence>0
	  } // end if ToCell != NULL
	} // end for i_col
  } // end for i_row

  // Free the allocated array
  FreeDouble2DArray(d_rel_Z_Vp);
  FreeDouble1DArray(d_summed_prob);

} // end Implement_Adaptive_Parameter_Averaging
//------------------------------------------------------------------------------
// Now make all newly colonised grid cells actually occupied//////////////
void convert_colonised_cellls_to_occupied(TRaft* Raft,
										  int arg_n_rows,
										  int arg_n_cols)
{
  int i_row, i_col;
  short int s_new;
  short int s_one = 1;
  short int s_zero = 0;

  TGridCell_ptr FromCell;

  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  FromCell = Raft->GetCell(i_row,
							   i_col);
	  if(FromCell != NULL) // Check if this cell is land
		{
		// check if the species has just dispersed to this cell
		s_new = FromCell->GetNewlyColonised();
		if(s_new>0) // the species is newly colonised
		  {
		  FromCell->SetOccurrence(s_one);
		  FromCell->SetNewlyColonised(s_zero);
		  } // end if s_new>0
		} // end if ToCell != NULL
	  } // end for i_col
	} // end for i_row
} // end convert_colonised_cellls_to_occupied

//******************************************************************************
//* A new function to generate a random deviate between a specified max and min
//* Written by Sam Moskwa
//******************************************************************************
float uniformRandomDeviate(float arg_f_min,
							float arg_f_max)
{
  float f_rand_deviate;

  f_rand_deviate =  arg_f_min + (arg_f_max-arg_f_min) * (float)rand()/(float)RAND_MAX;

  return f_rand_deviate;

} // end uniformRandomDeviate

/**************************************************************************
 * Generate Circular Fill
 * Even fills circle or doughnut
 * Should be nice to use with a cut out bit in the middle (arg_i_min_radius)
 * which can be done separately?
 * ARGUMENTS: const int arg_i_max_radius,       how far out
 *            const int arg_i_min_radius,       excluded zone to make a doughnut
 *            short ** arg_si_rel_grid,         output grid
 *************************************************************************
 int GenerateCircularFill(const int arg_i_max_radius,
						  const int arg_i_min_radius,
						  short ** arg_si_rel_grid)
{
int i_x, i_y;
int j_x, j_y;
int i_width;
int i_radius;
double d_radius;
int i_count;

  i_width= 2*arg_i_max_radius+2;
  i_radius=arg_i_max_radius+1;
  j_x=i_radius+i_x;
  j_y=i_radius+i_y;
  i_count=0;
  for(i_x=-i_radius;i_x<=i_radius;i_x++)
    {
    for(i_y=-i_radius;i_y<=i_radius;i_y++)
      {
	  d_radius=sqrt((double) ((i_x*i_x)+(i_y*i_y)));
      if(d_radius>=arg_i_min_radius&& d_radius<=arg_i_max_radius)
        {
        j_x=i_radius+i_x;
		j_y=i_radius+i_y;
        if(j_x>=0&&j_y>=0&&j_x<i_width&&j_y<i_width)
          {
          arg_si_rel_grid[j_x][j_y]=1;
          i_count++;
          }
        else
          j_y=1;
        }  // end if a valid point within radii
      } //  end for i_y
    } //  end for i_x
  return (i_count);
} // end func Generate Circular Fill
*/

//------------------------------------------------------------------------------
// 30 June 2015
// A function to calculate the Z & Vp (attribute mean and it's variance) for
// the newly colonised grid cell, based on those grid cells contributing
// propagules.
// A NEW FUNCTION BASED ON VARIANCE MIXTURE
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
											double** arg_d_rel_Vp)
{
  int i_env, i_rel;
  int i_count;
  double d_site_A_prob, d_site_B_prob;
  double d_site_A_Z, d_site_B_Z;
  double d_site_A_Vp, d_site_B_Vp;
  double d_A_proportion, d_B_proportion;
  double d_combined_Z;
  double d_combined_Vp;

  // Pointers to GridCells
  TGridCell_ptr ToCell;

  //pointers to local arrays
  double* d_Z_env_array;
  double* d_Vp_env_array;

  // Get the focal cell
  ToCell = Raft->GetCell(arg_i_row,
						 arg_i_col);
  if(ToCell != NULL) // Double check this cell is land
	{
	d_Z_env_array = ToCell->GetCurrentZ();
	d_Vp_env_array = ToCell->GetCurrentVp();
	// Loop through each env for the focal cell. If it's an adaptive env, calculate
	// the mean (z) & variance (Vp)
	for(i_env=0;i_env<arg_n_environments;i_env++)
	  {
	  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
		{
		// initialise the catching data
		i_count = 0;
		for(i_rel=0; i_rel<arg_n_rel; i_rel++)
		  {
		  if(arg_d_rel_Z[i_rel][i_env] > -9999)
			{
			if(i_count<1) // first occupied rel cell
			  {
			  d_site_A_prob = arg_d_rel_prob[i_rel];
			  d_site_A_Z = arg_d_rel_Z[i_rel][i_env];
			  d_site_A_Vp = arg_d_rel_Vp[i_rel][i_env];
			  //if(d_site_A_prob>0)   //NEW June 2016
			  //	i_count += 1;       //NEW June 2016
			  } // end if i_count<1
			else // other occupied rel cells
			  {
			  // catch the info for site B
			  d_site_B_prob = arg_d_rel_prob[i_rel];
			  d_site_B_Z = arg_d_rel_Z[i_rel][i_env];
			  d_site_B_Vp = arg_d_rel_Vp[i_rel][i_env];
			  // calculate the relative probability of site A & B
			  d_A_proportion = d_site_A_prob / (d_site_A_prob + d_site_B_prob);
			  d_B_proportion = d_site_B_prob / (d_site_A_prob + d_site_B_prob);
			  // calculate the combined mean
			  d_combined_Z = (d_site_A_Z * d_A_proportion) + (d_site_B_Z * d_B_proportion);
			  // and calculate the new combined (mixture) variance
			  d_combined_Vp = (d_A_proportion * d_site_A_Vp)+(d_B_proportion * d_site_B_Vp) +
							  ((d_A_proportion * (d_site_A_Z * d_site_A_Z)) +
							   (d_B_proportion * (d_site_B_Z * d_site_B_Z)) -
							   (((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z)) *
							   ((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z))));
			  // Now write the new mixture (combined) mean and variance to site A for the next
			  // increment through the relative cells
			  d_site_A_prob = d_site_A_prob + d_site_B_prob;
			  d_site_A_Z = d_combined_Z;
			  d_site_A_Vp = d_combined_Vp;
			  } // end else i_count<1
			i_count += 1;
			} // end if arg_d_rel_Z[i_rel][i_env] > -9999
		  } // end for i_rel
		// ensure the new Z doesn't exceed the evolutionary limit
		if(arg_b_env_low_adaptation[i_env]>0)
		  {
		  if(d_site_A_Z<arg_d_adapt_limit[i_env])
			d_site_A_Z = arg_d_adapt_limit[i_env];
		  if(d_site_A_Z > arg_f_env_low_fund_limit[i_env])   // A DODGY FIX - THIS WILL DESTROY GRID-BASED INITIAL Z
			d_site_A_Z = arg_f_env_low_fund_limit[i_env];
		  } // end if arg_b_env_low_adaptation[i_env]>0
		else
		  {
		  if(d_site_A_Z>arg_d_adapt_limit[i_env])
			d_site_A_Z = arg_d_adapt_limit[i_env];
		  if(d_site_A_Z < arg_f_env_high_fund_limit[i_env])   // A DODGY FIX - THIS WILL DESTROY GRID-BASED INITIAL Z
			d_site_A_Z = arg_f_env_high_fund_limit[i_env];
		  } // end else arg_b_env_low_adaptation[i_env]>0
		// Put the calculated mean and variance in the right spot
		d_Z_env_array[i_env] = d_site_A_Z;
		d_Vp_env_array[i_env] = d_site_A_Vp;
		} // end if arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0
	  } // end for i_env
	} // end if ToCell != NULL

} // end Calculate_Newly_Colonised_Z_Vp_Mixture

//------------------------------------------------------------------------------
// 1 July 2015
// Apply adaptive parameter averaging (admixture)
// Accounts for dispersal between occupied grid cells resulting in averaging of
// mean attributes and variance of attribute
// A NEW FUNCTION BASED ON VARIANCE MIXTURE
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
													double* arg_d_adapt_limit)
{
  // declare local variables
  int i_row, i_col, i_rel, i_env;
  int i_relrow, i_relcol;
  int i_count;
  short int s_occurrence;
  short int s_rel_occurrence;
  double d_site_A_prob, d_site_B_prob;
  double d_site_A_Z, d_site_B_Z;
  double d_site_A_Vp, d_site_B_Vp;
  double d_A_proportion, d_B_proportion;
  double d_combined_Z;
  double d_combined_Vp;

  TGridCell_ptr ToCell;
  TGridCell_ptr FromCell;

  // pointers to arrays
  double* d_Z_env_array;
  double* d_Vp_env_array;
  double* d_next_Z_env_array;
  double* d_next_Vp_env_array;
  double* d_rel_Z_env_array;
  double* d_rel_Vp_env_array;

// Loop through all the cells on the grid, calculate and make the correction
for(i_row=0; i_row<arg_n_rows; i_row++)
  {
  for(i_col=0; i_col<arg_n_cols; i_col++)
	{
	ToCell = Raft->GetCell(i_row,
						   i_col);
	if(ToCell != NULL) // Check if this cell is land
	  {
	  // check if the species occurs in this cell
	  s_occurrence = ToCell->GetOccurrence();
	  if(s_occurrence>0) // the species is present
		{
		// obtain the arrays for this grid cell
		d_Z_env_array = ToCell->GetCurrentZ();
		d_Vp_env_array = ToCell->GetCurrentVp();
		d_next_Z_env_array = ToCell->GetNextZ();
		d_next_Vp_env_array = ToCell->GetNextVp();
		for(i_env=0;i_env<arg_n_environments;i_env++)
		  {
		  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			{
			i_count=0;
			for(i_rel=0; i_rel<arg_n_rel; i_rel++)
			  {
			  i_relrow = i_row + arg_i_rel_row_col[i_rel][0];
			  i_relcol = i_col + arg_i_rel_row_col[i_rel][1];
			  // ensure the relative cell is still on the raft
			  if(i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols)
				{
				// check if this is in the domain
				FromCell = Raft->GetCell(i_relrow,
										 i_relcol);
				if(FromCell != NULL) // Check if this cell is land
				  {
				  // check if the species occurs in this cell
				  s_rel_occurrence = FromCell->GetOccurrence();
				  if(s_rel_occurrence>0) // the species is present
					{
					d_rel_Z_env_array = FromCell->GetCurrentZ();
					d_rel_Vp_env_array = FromCell->GetCurrentVp();
					if(i_count<1) // first occupied rel cell
					  {
					  d_site_A_prob = arg_d_rel_prob[i_rel];
					  d_site_A_Z = d_rel_Z_env_array[i_env];
					  d_site_A_Vp = d_rel_Vp_env_array[i_env];
					  if(d_site_A_prob>0)
						i_count += 1;
					  } // end if i_count<1
					else // other occupied rel cells
					  {
					  // catch the info for site B
					  d_site_B_prob = arg_d_rel_prob[i_rel];
					  d_site_B_Z = d_rel_Z_env_array[i_env];;
					  d_site_B_Vp = d_rel_Vp_env_array[i_env];
					  if(d_site_B_prob>0)
						{
						// calculate the relative probability of site A & B
						d_A_proportion = d_site_A_prob / (d_site_A_prob + d_site_B_prob);
						d_B_proportion = d_site_B_prob / (d_site_A_prob + d_site_B_prob);
						// calculate the combined mean
						d_combined_Z = (d_site_A_Z * d_A_proportion) + (d_site_B_Z * d_B_proportion);
						// and calculate the new combined (mixture) variance
						d_combined_Vp = (d_A_proportion * d_site_A_Vp)+(d_B_proportion * d_site_B_Vp) +
									  ((d_A_proportion * (d_site_A_Z * d_site_A_Z)) +
									   (d_B_proportion * (d_site_B_Z * d_site_B_Z)) -
									   (((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z)) *
									   ((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z))));
						// Now write the new mixture (combined) mean and variance to site A for the next
						// increment through the relative cells
						d_site_A_prob = d_site_A_prob + d_site_B_prob;
						d_site_A_Z = d_combined_Z;
						d_site_A_Vp = d_combined_Vp;
						i_count += 1;
						} // end if d_site_B_prob>0
					  } // end else i_count<1
					} // end if s_rel_occurrence>0
				  } // end if FromCell != NULL
				} // end if i_relrow>-1 && i_relrow<arg_n_rows && i_relcol>-1 && i_relcol<arg_n_cols
			  } // end for i_rel
			// Now mix the neighbouring cells and the focal cell, accounting for
			// the greater weighting the focal cell has in having a larger resident
			// population as opposed to a few successful dispersing propagules
			if(i_count > 0)
			  {
			  // catch the info for site B
			  d_site_B_prob = arg_d_resident_weighting;
			  d_site_B_Z = d_Z_env_array[i_env];
			  d_site_B_Vp = d_Vp_env_array[i_env];
			  // calculate the relative probability of site A & B
			  d_A_proportion = d_site_A_prob / (d_site_A_prob + d_site_B_prob);
			  d_B_proportion = d_site_B_prob / (d_site_A_prob + d_site_B_prob);
			  // calculate the combined mean
			  d_combined_Z = (d_site_A_Z * d_A_proportion) + (d_site_B_Z * d_B_proportion);
			  // and calculate the new combined (mixture) variance
			  d_combined_Vp = (d_A_proportion * d_site_A_Vp)+(d_B_proportion * d_site_B_Vp) +
							  ((d_A_proportion * (d_site_A_Z * d_site_A_Z)) +
							   (d_B_proportion * (d_site_B_Z * d_site_B_Z)) -
							   (((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z)) *
							   ((d_A_proportion * d_site_A_Z) + (d_B_proportion * d_site_B_Z))));
			  // and write the admixture derived Z & Vp values to the grid for the
			  // next time step, ensuring the new Z doesn't exceed the evolutionary limit
			  if(arg_b_env_low_adaptation[i_env]>0)
				{
				if(d_combined_Z<arg_d_adapt_limit[i_env])
				  d_combined_Z = arg_d_adapt_limit[i_env];
				} // end if arg_b_env_low_adaptation[i_env]>0
			  else
				{
				if(d_combined_Z>arg_d_adapt_limit[i_env])
				  d_combined_Z = arg_d_adapt_limit[i_env];
				} // end else arg_b_env_low_adaptation[i_env]>0
			  // & record the result
			  d_next_Z_env_array[i_env] = d_combined_Z;
			  d_next_Vp_env_array[i_env] = d_combined_Vp;
//TEMP
if(d_combined_Vp>1.106)
  i_count=i_count;
// END TEMP
			  } // end if i_count > 0
			else
			  {
			  d_next_Z_env_array[i_env] = d_Z_env_array[i_env];
			  d_next_Vp_env_array[i_env] = d_Vp_env_array[i_env];
			  } // end else i_count > 0
			} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
		  } // end for i_env
		} // end if s_occurrence>0
	  } // end if ToCell != NULL
	} // end for i_col
  } // end for i_row

// Now zip through the grid again and put the new Z & Vp values in the right spot
for(i_row=0; i_row<arg_n_rows; i_row++)
  {
  for(i_col=0; i_col<arg_n_cols; i_col++)
	{
	ToCell = Raft->GetCell(i_row,
						   i_col);
	if(ToCell != NULL) // Check if this cell is land
	  {
	  // check if the species occurs in this cell
	  s_occurrence = ToCell->GetOccurrence();
	  if(s_occurrence>0) // the species is present
		{
		// obtain the arrays for this grid cell
		d_Z_env_array = ToCell->GetCurrentZ();
		d_Vp_env_array = ToCell->GetCurrentVp();
		d_next_Z_env_array = ToCell->GetNextZ();
		d_next_Vp_env_array = ToCell->GetNextVp();
		for(i_env=0;i_env<arg_n_environments;i_env++)
		  {
		  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			{
			d_Z_env_array[i_env] = d_next_Z_env_array[i_env];
			d_Vp_env_array[i_env] = d_next_Vp_env_array[i_env];
			} // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
		  } // end for i_env
		} // end if s_occurrence>0
	  } // end if ToCell != NULL
	} // end for i_col
  } // end for i_row

} // end Implement_Adaptive_Parameter_Averaging
//------------------------------------------------------------------------------
// 7 July 2015
// A new function to summarize the distribution of this species at the end of
// the time point
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
								  int* arg_i_n_initial_occupied_ptr)
{
  // local variables
  int i_row, i_col, i_env;
  int i_occupied;
  int i_initial_occupied;
  short s_occurrence;

  TGridCell_ptr ToCell;

  // pointers to arrays
  double* d_Z_env_array;
  double* d_Vp_env_array;

  // initialise the catching arrays
  for(i_env=0;i_env<arg_n_environments;i_env++)
	{
	arg_d_mean_Z_env[i_env] = 0;
	arg_d_mean_Vp_env[i_env] = 0;
	} // end for i_env

  i_occupied = 0;
  i_initial_occupied = 0;
  for(i_row=0; i_row<arg_n_rows; i_row++)
	{
	for(i_col=0; i_col<arg_n_cols; i_col++)
	  {
	  ToCell = Raft->GetCell(i_row,
							   i_col);
	  if(ToCell != NULL) // Check if this cell is land
		{
		// check if the species occurs in this cell
		s_occurrence = ToCell->GetOccurrence();
		if(s_occurrence>0) // the species is present
		  {
		  i_occupied += 1;
		  // was the species in this cell at the beginning?
		  if(arg_s_initial_occurrence[i_row][i_col]>0)
			i_initial_occupied += 1;
		  // obtain the arrays for this grid cell
		  d_Z_env_array = ToCell->GetCurrentZ();
		  d_Vp_env_array = ToCell->GetCurrentVp();
		  for(i_env=0;i_env<arg_n_environments;i_env++)
			{
			if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
			  {
			  arg_d_mean_Z_env[i_env] += d_Z_env_array[i_env];
			  arg_d_mean_Vp_env[i_env] += d_Vp_env_array[i_env];
			  } // end if b_env_low_adaptation[i_env]>0 || b_env_high_adaptation[i_env]>0
			} // end for i_env
		  } // end if s_occurrence>0
		} // end if ToCell != NULL
	  } // end for i_col
	} // end for i_row

  // calculate the Z & Vp means
  if(i_occupied>0)
	{
	for(i_env=0;i_env<arg_n_environments;i_env++)
	  {
	  if(arg_b_env_low_adaptation[i_env]>0 || arg_b_env_high_adaptation[i_env]>0)
		{
		arg_d_mean_Z_env[i_env] = arg_d_mean_Z_env[i_env] / (double) (i_occupied);
		arg_d_mean_Vp_env[i_env] = arg_d_mean_Vp_env[i_env] / (double) (i_occupied);
		}// end if
	  } // end for i_env
	} // end if i_occupied>0

  // and write the number of cells occupied back to the arguments
  *arg_i_n_occupied_ptr = i_occupied;
  *arg_i_n_initial_occupied_ptr = i_initial_occupied;

} // end Generate_occurrence_summary
//------------------------------------------------------------------------------
