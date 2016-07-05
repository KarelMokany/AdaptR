//---------------------------------------------------------------------------
//*************************************************************
// MURU COMPACTOR
// Set of routines to compact layers of transformed gdm grids
// with a form of run length encoding  and shorts
// And hopefully some routines to get the data back out
//

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "DynamicArray.h"
#include "CompactorFunctions.h"

using namespace std;

//************************************************************
// Convert Grid Stack To Compact FLOAT Binary File
// takes a lot of binary grid files, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
int ConvertASCIIGridStackToCompactFloatBinaryFile(int i_nothing,
												  string arg_str_parameter_filename)
{
// Allocate local variables
int i_layer, i_step; // layer & step counter
int i_result_code=1; // good unless we hear otherwise!
int n_layers;
int i_pos; // position we have reached
int i_short_size;
int i_length;
int i_float_size;
int i_data_format;
int i_check;
int n_steps; // number of environment steps to loop through
int i_value;
int i_buf;
int i_do_nothing;
long int n_rows,n_cols; //local copies of arguments
long int i_row, i_col;  // row and column counters
long int i_start_run, i_end_run;
short i_s;
short int si_run_length;
unsigned short i_us;
float f_value;
double d_spos; // total offset of the start of each new line.
double d_value;

// Pointers
short* i_s_ptr;
unsigned short* i_us_ptr;
int* i_ptr;
long int* l_ptr;
float* f_ptr;    // for reading data in
float** f_data;
double* d_ptr;

// File streams
ifstream ProjectFile;
ifstream Out_namefile_stream;
ifstream* InFileStream_ptr;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;

// strings
string s_path,s_index_path;
string s_buffer;
string str_input_dir;
string str_output_dir;
string str_root_name_file;
string str_filename;
string str_env_name_filepath;
string str_ouput_namefile_path;
string str_ouput_file_path;
string str_out_path_base;

  // Get system data sizes
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);

// String arrays to hold file names
string str_layer_name_files[50]; // must be less than 50 env layers at each time step
string str_all_layer_filenames[5000][50]; // must be less than 5000 time steps
string str_output_grid_names[5000];

  // Read in the parameter file ////////////////////////////////////////////////
  ProjectFile.open(arg_str_parameter_filename.c_str());
  if(ProjectFile.is_open())
	{
    ProjectFile>>s_buffer;   // NCOLS
	ProjectFile>>n_cols;
	ProjectFile>>s_buffer;   // NROWS
	ProjectFile>>n_rows;
	ProjectFile>>s_buffer;   // NUM_ENV_LAYERS
	ProjectFile>>n_layers;
	ProjectFile>>s_buffer;   // NUM_ENV_STEPS
	ProjectFile>>n_steps;
	ProjectFile>>s_buffer;   // ENV_NAME_FILES
	i_buf=ProjectFile.get(); // read blank space
        for(i_layer=0; i_layer<n_layers; i_layer++)
	  {
	  getline(ProjectFile,s_buffer);
	  str_layer_name_files[i_layer]=s_buffer.c_str(); // read rest of line 
	  } // end for i_layer
//OLD	ProjectFile>>s_buffer;	// INPUT_FOLDER
//OLD	ProjectFile>>str_input_dir;
//OLD	ProjectFile>>s_buffer;	// OUTPUT_FOLDER
//OLD	ProjectFile>>str_output_dir;
	ProjectFile>>s_buffer;	// OUTPUT_NAME_FILE
	i_buf=ProjectFile.get(); // read blank space
	getline(ProjectFile,s_buffer);
	str_root_name_file=s_buffer.c_str(); // read rest of line;	
	ProjectFile>>s_buffer;  // PERCENTAGE_SUCCESS ";	
	ProjectFile>>i_check; 	// should be 100!
	ProjectFile.close();
	if(i_check!=100)
	  return(-2);
	} // end if ProjectFile.open
  else
	return(-1);
  //////////////////////////////////////////////////////////////////////////////



 
  // Now read in the environment grid input filenames //////////////////////////
  for(i_layer=0; i_layer<n_layers; i_layer++)
	{
	// read & open the name file
	s_buffer=str_layer_name_files[i_layer];
	str_filename=s_buffer.c_str();
//OLD	str_env_name_filepath = str_input_dir + str_filename;
	str_env_name_filepath = str_filename; //NEW
	ifstream Env_Name_FileStream;
	Env_Name_FileStream.open(str_env_name_filepath.c_str());
	for(i_step=0; i_step<n_steps; i_step++)
	  {
	  getline(Env_Name_FileStream,s_buffer);
	  str_all_layer_filenames[i_step][i_layer]=s_buffer.c_str(); // read rest of line	  
	  } // end for i_step
	Env_Name_FileStream.close();
	} // end for i_layer
  //////////////////////////////////////////////////////////////////////////////
 
  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat2DArray(n_layers,
							   n_cols);

  //----------------------------------------------------------------------------
//TEMP
/*ofstream Test_out_FileStream;//TEMP
string f_name; //TEMP
f_name = string("C:/Users/mok010/Karels Files/Mokany Files/My R applications/AdaptR/source/test.txt");//TEMP
Test_out_FileStream.open(f_name.c_str());//TEMP
Test_out_FileStream << "Start test file" << endl;//TEMP
Test_out_FileStream << "Output namefiles" << endl;//TEMP
str_ouput_namefile_path=str_root_name_file;//NEW
Test_out_FileStream << str_ouput_namefile_path << endl;//TEMP
Out_namefile_stream.open(str_ouput_namefile_path.c_str());
Test_out_FileStream << str_filename << endl;//TEMP	
Out_namefile_stream>>s_buffer;
Test_out_FileStream << s_buffer << endl;//TEMP
str_ouput_file_path = s_buffer+".mcg";//NEW  
Test_out_FileStream << str_ouput_file_path << endl;//TEMP
Out_namefile_stream>>s_buffer;//TEMP
Test_out_FileStream << s_buffer << endl;//TEMP
str_ouput_file_path = s_buffer+".mcg";//NEW  
Test_out_FileStream << str_ouput_file_path << endl;//TEMP
s_buffer=str_all_layer_filenames[0][0].c_str();//TEMP
str_filename=s_buffer;//TEMP
Test_out_FileStream << str_filename << endl; //TEMP
Test_out_FileStream << "End test file" << endl;//TEMP
Test_out_FileStream.close();//TEMP
return (33);//TEMP
*/
// END TEMP
  
  // Now loop through all the time steps, compact the environment grids & write
  // them out to binary
  // First open the name file for the output names
//OLD  str_ouput_namefile_path=str_input_dir+str_root_name_file+".txt";
  str_ouput_namefile_path=str_root_name_file;//NEW	 
  Out_namefile_stream.open(str_ouput_namefile_path.c_str());
  for(i_step=0; i_step<n_steps; i_step++)
	{
    // Open up all grid files to compact
	InFileStream_ptr = new ifstream[n_layers];
	for(i_layer=0;i_layer<n_layers;i_layer++)
	  {
	  s_buffer=str_all_layer_filenames[i_step][i_layer].c_str();
	  str_filename=s_buffer; 
	  InFileStream_ptr[i_layer].open(str_filename.c_str());
// NEW // check opening of input files //????????????????????????????????????????????????//
/*
	  if(InFileStream_ptr[i_layer].is_open())
	    { 
	     i_do_nothing = 1;
	    } // end if InFileStream_ptr[i_layer].is_open()
	  else
	    {
		delete [] InFileStream_ptr;
		FreeFloat2DArray(f_data);
		return (-3);	
	    } // end else InFileStream_ptr[i_layer].is_open()
*/
// END NEW // check opening of input files //??????????????????????????????????????????//	
	  // read in the header
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>i_value;
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>i_value;
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>d_value;
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>d_value;
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>d_value;
	  InFileStream_ptr[i_layer]>>s_buffer;
	  InFileStream_ptr[i_layer]>>d_value;  // changed from i_value, incase we get float no data
	  } // end for i_layer
	  
	// Open big results file AS BINARY
//OLD	Out_namefile_stream>>s_buffer;
//OLD	str_ouput_file_path = str_output_dir+s_buffer+".mcg";
	getline(Out_namefile_stream,s_buffer);
    str_out_path_base = s_buffer.c_str();
	str_ouput_file_path = str_out_path_base + ".mcg";//NEW
	ofstream OutFileStream;
	OutFileStream.open(str_ouput_file_path.c_str(), ios::out |ios::binary);
// NEW // check opening of output files //??????????????????????????????????????????????//
/*
	  if(OutFileStream.is_open())
	    { 
	     i_do_nothing = 1;
	    } // end if OutFileStream.is_open()
	  else
	    {
		delete [] InFileStream_ptr;
		FreeFloat2DArray(f_data);
		return (-4);
		} // end else OutFileStream.is_open()
*/
// END NEW // check opening of output files //???????????????????????????????????????????//	
	// build index file name ( to hold seekposes for the beginning of each line)
	//s_path=str_ouput_file_path; //OLD
	//i_length=s_path.length(); //OLD
	//s_path.resize(i_length-4); //OLD
	//s_index_path=s_path+".mci"; //OLD
	s_index_path = str_out_path_base + ".mci";
	ofstream IndexFileStream;
	IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);
// NEW // check opening of output files //???????????????????????????????????????????????//
/*
	  if(IndexFileStream.is_open())
	    { 
	     i_do_nothing = 1;
	    } // end if IndexFileStream.is_open()
	  else
	    {
		delete [] InFileStream_ptr;
		FreeFloat2DArray(f_data);
		return (-5);
		} // end else IndexFileStream.is_open()
*/
// END NEW // check opening of output files //?????????????????????????????????.????????//
	// build header file name ( to hold ncols, nrows, nlayers )
//OLD	s_path=str_output_dir+s_buffer+".mch";
	//s_path=str_ouput_file_path+".mch"; //NEW
	s_path = str_out_path_base + ".mch"; //NEW
	ofstream HeaderFileStream;
	HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
// NEW // check opening of output files //???????????????????????????????????????????????//
/*
	  if(HeaderFileStream.is_open())
	    { 
	     i_do_nothing = 1;
	    } // end if IndexFileStream.is_open()
	  else
	    {
		delete [] InFileStream_ptr;
		FreeFloat2DArray(f_data);
		return (-6);
		} // end else IndexFileStream.is_open()
*/
// END NEW // check opening of output files //?????????????????????????????????.????????//

	l_ptr=&n_cols;
	HeaderFileStream.write((char*)l_ptr,sizeof(long int));
	l_ptr=&n_rows;
	HeaderFileStream.write((char*)l_ptr,sizeof(long int));
	i_ptr=&n_layers;
	HeaderFileStream.write((char*)i_ptr,sizeof(int));
	i_data_format=1;
	i_ptr=&i_data_format;
	HeaderFileStream.write((char*)i_ptr,sizeof(int));
	HeaderFileStream.close();

	d_spos=0; // zero seekpos counter
	//Read whole lines into buffers
	for(i_row=0;i_row<n_rows;i_row++)
	  {
	  // read in whole line and process
	  for(i_col=0;i_col<n_cols;i_col++)
		{
		// Read in values for each cell
		for(i_layer=0;i_layer<n_layers;i_layer++)
		  {
		  InFileStream_ptr[i_layer]>>d_value;
                   f_data[i_layer][i_col]= (float) d_value;
                  //InFileStream_ptr[i_layer]>>f_data[i_layer][i_col];  // commented out for checking
		  }// end for i_layer
		} // end for i_col
	  // Assume all layers are aligned, so we need only do fiddling in the topmost layer
	  i_pos=0;
	  i_col=0;
	  //Work through the whole row
	  while(i_pos<n_cols)
		{
		// split into nodata or data runs
		if(f_data[0][i_pos] < -9000)
		  {
		  i_start_run=i_pos;
		  si_run_length=0;
		  while(f_data[0][i_pos] < -9000 && i_pos < n_cols)
			{
			si_run_length++;  // 1 first time through
			i_pos++;
			if(si_run_length >= 32700) // size of short: don't run off the end of a short for really big grids
			  break;
			}// end while still nodata
		  // we have found the end of the run.
		  // Since this is no_data.. convert to negative
		  si_run_length=0-si_run_length;
		  // Write out this only
		  i_s_ptr=&(si_run_length);
		  OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
		  d_spos+=i_short_size;
		  } // end if NODATA run
		else
		  {
		  i_start_run=i_pos;
		  si_run_length=0;
		  while(f_data[0][i_pos] >= -9000 && i_pos<n_cols)
			{
			si_run_length++;
			i_pos++;
			if(si_run_length >= 32700) // size of short
			  break;
			}// end while still data
		  // we have found the end of the run.

		  //Go back to start of run and read out ALL LAYERS
		  // Write out run length first
		  i_s_ptr=&(si_run_length);
		  OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
		  d_spos+=i_short_size;
		  // Read out whoel run, all layer data for each point
		  // Output data for each Grid cell for all layers
		  // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
		  // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
		  i_end_run=i_start_run+(int)si_run_length;
		  for(i_col=i_start_run;i_col<i_end_run;i_col++)
			{
			for(i_layer=0;i_layer<n_layers;i_layer++)
			  {
			  f_value=f_data[i_layer][i_col];
			  f_ptr=&(f_value);
			  OutFileStream.write((char*)f_ptr,i_float_size);   // unshortened value for this layer
			  d_spos+=i_float_size;
			  }// end for i_layer
			}// end for i_col start to end run
		  } // end else DATA run
		} // end while(i_pos<n_cols)
	  // At the end of the row, output the offset from the start of the file
	  d_ptr=&d_spos;
	  IndexFileStream.write((char*)d_ptr,sizeof(double));

	  }// end for i_row
	// Output an end of file token
	i_s = -32000;
	i_s_ptr=&i_s;
	OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
	//Close files
	OutFileStream.close();
	IndexFileStream.close();
	for(i_layer=0;i_layer<n_layers;i_layer++)
	  InFileStream_ptr[i_layer].close();

	// delete filestream pointer
	delete [] InFileStream_ptr;
	} // end for i_step

  // Free up arrays
  FreeFloat2DArray(f_data);

  return i_result_code;
}// end ConvertGridStackToCompactFLOATBinaryFile


/*

//************************************************************
// Convert Condition Grid  To Compact FLOAT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
unsigned short* i_us_ptr;
unsigned short i_us;
float f_value;
string s_path,s_index_path;
int i_short_size, i_unshort_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_FLOAT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
      f_ptr=&(f_data[i_col]);
      //Read to f_data array
      InFileStream.read((char*)f_ptr,sizeof(float));
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;  // 1 first time through
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
        si_run_length=0-si_run_length;
        // Write out this only
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          f_ptr=&f_value;
          OutFileStream.write((char*)f_ptr,i_float_size);   // value for this layer
          d_spos+=i_float_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  f_value=K_EOF;
  f_ptr=&f_value;
  OutFileStream.write((char*)f_ptr,i_float_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertConditionGridToCompactBinaryFile

//************************************************************
// Convert an ASCII Condition Grid  To Compact FLOAT Binary File
// takes a ascii grid file (*.asc), and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
extern int ConvertASCIIConditionGridToCompactFloatBinaryFile(string arg_s_condition_file,
                                                             string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
unsigned short* i_us_ptr;
unsigned short i_us;
float f_value;
string s_path,s_index_path;
int i_short_size, i_unshort_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_float_size;
int i_data_format;
double* d_ptr;
long int* l_ptr;
float f_indata;
  int i_nodata;
  double d_xllcorner, d_yllcorner;
  double d_cellsize;
  string str_buffer;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str());
  InFileStream>>str_buffer;
  InFileStream>>n_cols;
  InFileStream>>str_buffer;
  InFileStream>>n_rows;
  InFileStream>>str_buffer;
  InFileStream>>d_xllcorner;
  InFileStream>>str_buffer;
  InFileStream>>d_yllcorner;
  InFileStream>>str_buffer;
  InFileStream>>d_cellsize;
  InFileStream>>str_buffer;
  InFileStream>>i_nodata;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols);

  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_FLOAT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
	// read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
	  InFileStream >> f_indata;
	  f_data[i_col] = f_indata;
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
		si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
		  {
          si_run_length++;  // 1 first time through
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
        si_run_length=0-si_run_length;
        // Write out this only
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          f_ptr=&f_value;
          OutFileStream.write((char*)f_ptr,i_float_size);   // value for this layer
          d_spos+=i_float_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  f_value=K_EOF;
  f_ptr=&f_value;
  OutFileStream.write((char*)f_ptr,i_float_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertASCIIConditionGridToCompactBinaryFile
*/

//************************************************************
// Convert Richness Grid  To Compact FLOAT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are floats
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
//****************************************************************
int ConvertASCIIRichnessGridToCompactFloatBinaryFile(int i_nothing,
													  string arg_str_parameter_filename)
{

//Allocate local variables
int i_layer, i_step; // layer & step counter
int i_result_code=1; // good unless we hear otherwise!
int n_layers;
int i_pos; // position we have reached
int i_short_size;
int i_length;
int i_float_size;
int i_data_format;
int i_check;
int n_steps; // number of environment steps to loop through
int i_value;
int i_buf;
long int n_rows,n_cols; //local copies of arguments
long int i_row, i_col;  // row and column counters
long int i_start_run, i_end_run;
short i_s;
short int si_run_length;
unsigned short i_us;
float f_value;
double d_spos; // total offset of the start of each new line.
double d_value;

// Pointers
short* i_s_ptr;
unsigned short* i_us_ptr;
int* i_ptr;
long int* l_ptr;
float* f_ptr;    // for reading data in
float* f_data;
double* d_ptr;

// File streams
ifstream ProjectFile;
ifstream Out_namefile_stream;

// strings
string s_path,s_index_path;
string s_buffer;
string str_input_dir;
string str_output_dir;
string str_root_name_file;
string str_filename;
string str_env_name_filepath;
string str_ouput_namefile_path;
string str_ouput_file_path;

  // Get system data sizes
  i_short_size= sizeof(short);
  i_float_size=sizeof(float);

// String arrays to hold file names
string str_layer_name_files[30]; // must be less than 30 env layers at each time step
string str_all_layer_filenames[800][30]; // must be less than 800 time steps
string str_output_grid_names[800];

  // Read in the parameter file ////////////////////////////////////////////////
  ProjectFile.open(arg_str_parameter_filename.c_str());
  if(ProjectFile.is_open())
	{
	ProjectFile>>s_buffer;   // NCOLS
	ProjectFile>>n_cols;
	ProjectFile>>s_buffer;   // NROWS
	ProjectFile>>n_rows;
	ProjectFile>>s_buffer;   // NUM_ENV_LAYERS
	ProjectFile>>n_layers;
	ProjectFile>>s_buffer;   // NUM_ENV_STEPS
	ProjectFile>>n_steps;
	ProjectFile>>s_buffer;   // ENV_NAME_FILES
	i_buf=ProjectFile.get(); // read blank space
        for(i_layer=0; i_layer<n_layers; i_layer++)
	  {
	  getline(ProjectFile,s_buffer);
	  str_layer_name_files[i_layer]=s_buffer.c_str(); // read rest of line
	  } // end for i_layer
	ProjectFile>>s_buffer;	// INPUT_FOLDER
	ProjectFile>>str_input_dir;
	ProjectFile>>s_buffer;	// OUTPUT_FOLDER
	ProjectFile>>str_output_dir;
	ProjectFile>>s_buffer;	// OUTPUT_ROOT_NAME_FILE
	ProjectFile>>str_root_name_file;
	ProjectFile>>s_buffer;  // PERCENTAGE_SUCCESS ";
	ProjectFile>>i_check; 	// should be 100!
	ProjectFile.close();
	if(i_check!=100)
	  return(0);
	} // end if ProjectFile.open
  else
	return(0);
  //////////////////////////////////////////////////////////////////////////////

  // Now read in the environment grid input filenames //////////////////////////
  // read & open the name file
  s_buffer=str_layer_name_files[0];
  str_filename=s_buffer.c_str();
  str_env_name_filepath = str_input_dir + str_filename;
  ifstream Env_Name_FileStream;
  Env_Name_FileStream.open(str_env_name_filepath.c_str());
  for(i_step=0; i_step<n_steps; i_step++)
	{
	getline(Env_Name_FileStream,s_buffer);
	str_all_layer_filenames[i_step][0]=s_buffer.c_str(); // read rest of line
	} // end for i_step
  Env_Name_FileStream.close();
  //////////////////////////////////////////////////////////////////////////////

  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat1DArray(n_cols);

  //----------------------------------------------------------------------------
  // Now loop through all the time steps, compact the environment grids & write
  // them out to binary
  // First open the name file for the output names
  str_ouput_namefile_path=str_input_dir+str_root_name_file+".txt";
  Out_namefile_stream.open(str_ouput_namefile_path.c_str());
  for(i_step=0; i_step<n_steps; i_step++)
	{
	// Open up the grid file to compact
	s_buffer=str_all_layer_filenames[i_step][0].c_str();
	str_filename=s_buffer;
	ifstream InFileStream;
	InFileStream.open(str_filename.c_str());

	// read in the header
	InFileStream>>s_buffer;
	InFileStream>>i_value;
	InFileStream>>s_buffer;
	InFileStream>>i_value;
	InFileStream>>s_buffer;
	InFileStream>>d_value;
	InFileStream>>s_buffer;
	InFileStream>>d_value;
	InFileStream>>s_buffer;
	InFileStream>>d_value;
	InFileStream>>s_buffer;
	InFileStream>>d_value;

	// Open results file AS BINARY
	Out_namefile_stream>>s_buffer;
	str_ouput_file_path = str_output_dir+s_buffer+".mcr";
	ofstream OutFileStream;
	OutFileStream.open(str_ouput_file_path.c_str(), ios::out |ios::binary);

	// build index file name ( to hold seekposes for the beginning of each line)
	s_path=str_ouput_file_path;
	i_length=s_path.length();
	s_path.resize(i_length-4);
	s_index_path=s_path+".mci";
	ofstream IndexFileStream;
	IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

	// build header file name ( to hold ncols, nrows, nlayers )
	s_path=str_output_dir+s_buffer+".mch";
	ofstream HeaderFileStream;
	HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
	l_ptr=&n_cols;
	HeaderFileStream.write((char*)l_ptr,sizeof(long int));
	l_ptr=&n_rows;
	HeaderFileStream.write((char*)l_ptr,sizeof(long int));
	i_ptr=&n_layers;
	n_layers=1;
	HeaderFileStream.write((char*)i_ptr,sizeof(int));
	i_data_format=1;
	i_ptr=&i_data_format;
	HeaderFileStream.write((char*)i_ptr,sizeof(int));
	HeaderFileStream.close();

	d_spos=0; // zero seekpos counter
	//Read whole lines into buffers
	for(i_row=0;i_row<n_rows;i_row++)
	  {
	  // read in whole line and process
	  for(i_col=0;i_col<n_cols;i_col++)
		{
		InFileStream>>d_value;
		f_data[i_col] = (float) (d_value);
		}
	  // Assume all layers are aligned, so we need only do fiddling in the topmost layer
	  i_pos=0;
	  i_col=0;
	  //Work through the whole row
	  while(i_pos<n_cols)
		{
		// split into nodata or data runs
		if(f_data[i_pos] < -9000)
		  {
		  i_start_run=i_pos;
		  si_run_length=0;
		  while(f_data[i_pos] < -9000 && i_pos < n_cols)
			{
			si_run_length++;  // 1 first time through
			i_pos++;
			if(si_run_length >= 32700) // size of short: don't run off the end of a short for really big grids
			  break;
			}// end while still nodata
		  // we have found the end of the run.
		  // Since this is no_data.. convert to negative
		  si_run_length=0-si_run_length;
		  // Write out this only
		  i_s_ptr=&(si_run_length);
		  OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
		  d_spos+=i_short_size;
		  } // end if NODATA run
		else
		  {
		  i_start_run=i_pos;
		  si_run_length=0;
		  while(f_data[i_pos] >= -9000 && i_pos<n_cols)
			{
			si_run_length++;
			i_pos++;
			if(si_run_length >= 32700) // size of short
			  break;
			}// end while still data
		  // we have found the end of the run.

		  //Go back to start of run and read out ALL LAYERS
		  // Write out run length first
		  i_s_ptr=&(si_run_length);
		  OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
		  d_spos+=i_short_size;
		  // Read out whole run, all layer data for each point
		  // Output data for each Grid cell for all layers
		  // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
		  // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
		  i_end_run=i_start_run+(int)si_run_length;
		  for(i_col=i_start_run;i_col<i_end_run;i_col++)
			{
			f_value=f_data[i_col];
			f_ptr=&(f_value);
			OutFileStream.write((char*)f_ptr,i_float_size);   // unshortened value for this layer
			d_spos+=i_float_size;
			}// end for i_col start to end run
		  } // end else DATA run
		} // end while(i_pos<n_cols)
	  // At the end of the row, output the offset from the start of the file
	  d_ptr=&d_spos;
	  IndexFileStream.write((char*)d_ptr,sizeof(double));
	  }// end for i_row

	// Output an end of file token
	i_s = -32000;
	i_s_ptr=&i_s;
	OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
	//Close files
	OutFileStream.close();
	IndexFileStream.close();
	InFileStream.close();
	} // end for i_step

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;

}// end ConvertASCIIRichnessGridToCompactFloatBinaryFile

/*

//************************************************************
// Convert Grid Stack To Compact Binary SHORT File
// takes a lot of binary grid files, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
// DOES NOT HANDLE NEGATIVE NUMBERS!
//****************************************************************
extern int ConvertGridStackToCompactShortBinaryFile(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
                                                    const int arg_i_num_layers,
                                                    string arg_s_layer_files[],
                                                    string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float** f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
long int i_pos; // position we have reached
ifstream* InFileStream_ptr;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
unsigned short* i_us_ptr;
unsigned short i_us;
float f_value;
string s_path,s_index_path;
int i_short_size, i_unshort_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
long int* l_ptr;
int i_buffer;
int i_data_format;
double* d_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_unshort_size= sizeof(unsigned short);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;
  n_layers= arg_i_num_layers;

  // Allocate a set of arrays to hold the lines for all layers
  f_data= AllocateFloat2DArray(n_layers,
                               n_cols);

  // Open up all grid files to compact
  InFileStream_ptr = new ifstream[n_layers];
  for(i_layer=0;i_layer<n_layers;i_layer++)
    {
    InFileStream_ptr[i_layer].open(arg_s_layer_files[i_layer].c_str(), ios::in |ios::binary);
    }
  // Open big results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_USHORT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in values for each cell
      for(i_layer=0;i_layer<n_layers;i_layer++)
        {
        f_ptr=&(f_data[i_layer][i_col]);
        //Read to f_data array
        InFileStream_ptr[i_layer].read((char*)f_ptr,sizeof(float));
        }// end for i_layer
      } // end for i_col

    // Assume all layers are aligned, so we need only
    // do fiddling in the topmost layer
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[0][i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
      //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[0][i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;  // 1 first time through
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
        si_run_length=0-si_run_length;
        // Write out this only
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[0][i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          for(i_layer=0;i_layer<n_layers;i_layer++)
            {
            f_value=f_data[i_layer][i_col];
            if(f_value> K_TRANS_MAX)//K_TRANS_MAX)
              {
              i_result_code=THREE;
              f_value= K_TRANS_MAX; // truncate to max value of 13.125
              }
            if(f_value<0)
              f_value=0;
            f_value*=K_SHORT_SCALE;
            i_us=(unsigned short)f_value;
            i_us_ptr=&(i_us);
            OutFileStream.write((char*)i_us_ptr,i_unshort_size);   // unshortened value for this layer
            d_spos+=i_unshort_size;
            }// end for i_layer
          }// end for i_col start to end run
        } // end DATA run
      }
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));

    }// end for i_row
  // Output an end of file token
  si_run_length=K_EOF;
  i_s_ptr=&(si_run_length);
  OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();

  for(i_layer=0;i_layer<n_layers;i_layer++)
     InFileStream_ptr[i_layer].close();

  // Free up arrays
  FreeFloat2DArray(f_data);
  delete [] InFileStream_ptr;

  return i_result_code;
}// end ConvertGridStackToCompactSHORTBinaryFile

//************************************************************
// Convert Condition Grid  To Compact SHORT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
// DOES NOT HANDLE NEGATIVE NUMBERS!
//****************************************************************
extern int ConvertConditionGridToCompactUShortBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
unsigned short* i_us_ptr;
unsigned short i_us;
float f_value;
string s_path,s_index_path;
int i_short_size, i_unshort_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_data_format;
double* d_ptr;
long int* l_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_unshort_size= sizeof(unsigned short);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_USHORT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
      f_ptr=&(f_data[i_col]);
      //Read to f_data array
      InFileStream.read((char*)f_ptr,sizeof(float));
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;  // 1 first time through
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
        si_run_length=0-si_run_length;
        // Write out this only
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_SHORT_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          if(f_value> K_TRANS_MAX)//K_TRANS_MAX)
            {
            i_result_code=THREE;
            f_value= K_TRANS_MAX; // truncate to max value of 13.125
            }
          if(f_value<0)
            f_value=0;
          f_value*=K_SHORT_SCALE;
          i_us=(unsigned short)f_value;
          i_us_ptr=&(i_us);
          OutFileStream.write((char*)i_us_ptr,i_unshort_size);   // unshortened value for this layer
          d_spos+=i_unshort_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  si_run_length=K_EOF;
  i_s_ptr=&(si_run_length);
  OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertConditionGridToCompactSHORTBinaryFile

//************************************************************
// Convert Richness Grid  To Compact SHORT Binary File
// takes a binary grid file in DIVA GIS .gri format, and squishes them a lot
// floats are converted to unsigned short ints
// the format is RUN LENGTH : datapoint 1 layer 1, datapoint1 layer 2 ..datapoint1 layer n:
// datapoint 2 layer 1, datapoint2 layer 2 ..datapoint2 layer n:  to datapoint n layer n
// NODATA runs are given as NEGATIVES, with no data afterwards.
// e.g. -5 5 34 23 24 25 21 23 25 12 25 25 -40 for a 2 layer grid.. is
// 5 nodata cells, 5 data cells, first cell, 34, 23 etc, followed by 40 nodatas
// DOES NOT HANDLE NEGATIVE NUMBERS!
//****************************************************************
extern int ConvertRichnessGridToCompactUShortBinaryFile( const long int arg_i_num_rows,
                                                         const long int arg_i_num_cols,
                                                         string arg_s_condition_file,
                                                         string arg_s_res_path)
{
int i_result_code;
long int n_rows,n_cols; //local copies of arguments
int n_layers;
long int i_row, i_col;  // row and column counters
int i_layer;       //layer/file counter
float* f_data;
//float** f_cell_pptr;     // Pointer to float pointer to make array of a pointer for each layer
int i_pos; // position we have reached
ifstream InFileStream;
ofstream OutFileStream;
ofstream IndexFileStream;
ofstream HeaderFileStream;
float* f_ptr;    // for reading data in
long int i_start_run, i_end_run;
short int si_run_length;
short* i_s_ptr;
short i_s;
unsigned short* i_us_ptr;
unsigned short i_us;
float f_value;
string s_path,s_index_path;
int i_short_size, i_unshort_size;
double d_spos; // total offset of the start of each new line.
int i_length;
int* i_ptr;
int i_data_format;
double* d_ptr;
long int* l_ptr;

  i_result_code=1; // good unless we hear otherwise!
  i_short_size= sizeof(short);
  i_unshort_size= sizeof(unsigned short);
  n_rows= arg_i_num_rows;
  n_cols= arg_i_num_cols;

  // Allocate an array to hold the lines
  f_data= AllocateFloat1DArray(n_cols);

  // Open up condition grid file to compact
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  // Open results file AS BINARY
  OutFileStream.open(arg_s_res_path.c_str(), ios::out |ios::binary);
  // build index file name ( to hold seekposes for the beginning of each line)
  s_path=arg_s_res_path;
  i_length=s_path.length();
  s_path.resize(i_length-4);
  s_index_path=s_path+".mci";
  IndexFileStream.open(s_index_path.c_str(), ios::out |ios::binary);

  // Quickly open and close and output header file details: cols. rows. layers.
  s_path=s_path+".mch";
  HeaderFileStream.open(s_path.c_str(), ios::out |ios::binary);
  l_ptr=&n_cols;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  l_ptr=&n_rows;
  HeaderFileStream.write((char*)l_ptr,sizeof(long int));
  i_ptr=&n_layers;
  n_layers=1;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  i_data_format=K_USHORT;
  i_ptr=&i_data_format;
  HeaderFileStream.write((char*)i_ptr,sizeof(int));
  HeaderFileStream.close();

  d_spos=0; // zero seekpos counter
  //Read whole lines into buffers
  for(i_row=0;i_row<n_rows;i_row++)
    {
    // read in whole line and process
    for(i_col=0;i_col<n_cols;i_col++)
      {
      // Read in value for each cell
      f_ptr=&(f_data[i_col]);
      //Read to f_data array
      InFileStream.read((char*)f_ptr,sizeof(float));
      } // end for i_col
    i_pos=0;
    i_col=0;
    //Work through the whole row
    while(i_pos<n_cols)
      {
      // split into nodata or data runs
      if(f_data[i_pos]<K_NO_DATA_ISH)
        {
        //NO_DATA
        //  i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]<K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;  // 1 first time through
          i_pos++;
          if(si_run_length>=K_RICH_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still nodata
        // we have found the end of the run.
        // Since this is no_data.. convert to negative
        si_run_length=0-si_run_length;
        // Write out this only
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // negative run length
        d_spos+=i_short_size;
        } // end NODATA run
      else
        {
        //DATA
       // i_col=i_pos;
        i_start_run=i_pos;
        si_run_length=0;
        while(f_data[i_pos]>=K_NO_DATA_ISH&&i_pos<n_cols)
          {
          si_run_length++;
          i_pos++;
          if(si_run_length>=K_RICH_MAX) // don't run off the end of a short for really big grids
            break;
          }// end while still data
        // we have found the end of the run.

        //Go back to start of run and read out ALL LAYERS
        // Write out run length first
        i_s_ptr=&(si_run_length);
        OutFileStream.write((char*)i_s_ptr,i_short_size);   // positive run length
        d_spos+=i_short_size;
        // Read out whoel run, all layer data for each point
        // Output data for each Grid cell for all layers
        // in order of the layers.. e.g. for 3 layers 1, 2,3, with run length 2
        // 2(runlength) 1 2 3 (firstcell) 1 2 3 (second cell)
        i_end_run=i_start_run+(int)si_run_length;
        for(i_col=i_start_run;i_col<i_end_run;i_col++)
          {
          f_value=f_data[i_col];
          if(f_value> K_RICH_MAX)//K_TRANS_MAX)
            {
            i_result_code=THREE;
            f_value= K_RICH_MAX; // truncate to max value of 327
            }
          if(f_value<0)
            f_value=0;
          f_value*=K_RICH_SCALE;
          i_us=(unsigned short)f_value;
          i_us_ptr=&(i_us);
          OutFileStream.write((char*)i_us_ptr,i_unshort_size);   // unshortened value for this layer
          d_spos+=i_unshort_size;
          }// end for i_col start to end run
        } // end DATA run
      } // end while
    // At the end of the row, output the offset from the start of the file
    d_ptr=&d_spos;
    IndexFileStream.write((char*)d_ptr,sizeof(double));
    }// end for i_row

  // Output an end of file token
  si_run_length=K_EOF;
  i_s_ptr=&(si_run_length);
  OutFileStream.write((char*)i_s_ptr,i_short_size);   // eoftoken
  //Close files
  OutFileStream.close();
  IndexFileStream.close();
  InFileStream.close();

  // Free up arrays
  FreeFloat1DArray(f_data);

  return i_result_code;
}// end ConvertRichnessGridToCompactFLOATBinaryFile

*/
//---------------------------------------------------------------------------

