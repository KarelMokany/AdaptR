//*****************************************************************************
// TRAFT . CPP:       (MURU)
// written by Tom Harwood Jan-Aug 2011
// Data structure to store data only for LAND, ( NULL for sea)
// and deal with multiple time steps and interpolation
// Designed to be nice and self contained so it has a realy simple interface
// Use of increment to interpolate between known data steps
// Can take allocated species data
// Also deals with dispersal or sampling over the grid
//
//*****************************************************************************
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "defines.h"
//#include "structures.h"
//#include "Generate.h"
#include "DynamicArray.h"
#include "GridCell.h"  // this contains a lot of the stuff we use here
#include "Raft.h"
using namespace std;
//****************************************************************************
// T RAFT: Constructor
// Creates an empty Raft, which is a grid nrows*ncols of TGridCell_ptr, all
// of which are set to NULL for sea/nodata for now
// We can then upload stuff to bits of it as required
//****************************************************************************
TRaft::TRaft(int arg_n_rows,
			 int arg_n_cols,
			 short** arg_s_condition)

{
  long int i_row; // row counter
  long int i_col;// col counter
  int i_land_count;
  //string s_land_file;
  //string s_sample_file;
  //string s_file;
  n_rows=arg_n_rows;
  n_cols=arg_n_cols;

  //Run_Parameters=arg_Run_Parameters;
  //n_rows=Run_Parameters->n_rows;
  //n_cols= Run_Parameters->n_cols;
  //n_layers=Run_Parameters->n_layers;
  // Step 1: Allocate the 2D array of TGridCellpointers
  // Allocate memory for array of elements of column ..
  // i.e. list of pointers to first element of each row
  Map = new TGridCell_ptr*[arg_n_rows];
  // Allocate memory for array of elements of all the rows
  // Do this as one big block to speed things up
  TGridCell_ptr* TCD_row_ptr = new TGridCell_ptr[arg_n_rows*arg_n_cols];
  // Now point the pointers in the right place
  for(i_row=0; i_row<arg_n_rows; i_row++)
    {
	*(Map+i_row)=TCD_row_ptr; // Could also be written as T_new_pptr[i_row]=
	TCD_row_ptr+=arg_n_cols;       // skip forward to the beginning of the next row
    }
  // Allocated.. Step 2.. Nullify
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  Map[i_row][i_col]=NULL;
      } // end for i_col
    } // end for i_row
  // Make space to hold the loaded row details
  n_land_in_row= AllocateLong1DArray(arg_n_rows);
  // Raft is now created
  // Create non NULL grid cells for land
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	i_land_count=0;
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(arg_s_condition[i_row][i_col]>0)
		{
		Map[i_row][i_col]=new TGridCell();
		// count valid cells
		i_land_count++;
		} // end if arg_s_condition[i_row][i_col]>0
	  } // end for i_col
	n_land_in_row[i_row]=i_land_count;
	} // end for i_row

 // initialise all the data time series and load into the raft
 //InitialiseAllSeries();

} // end constr 2

//---------------------------------------------------------------------------
   // destructor    // NEED TO REVISE THIS ********
TRaft::~TRaft(void)
{
long int i_row; // row counter
long int i_col;// col counter

  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        delete (Map[i_row][i_col]);
      }
    }
  // Free up both levels of the allocated array
  delete [] *Map;
  delete [] Map;

  FreeLong1DArray(n_land_in_row);

}

/*****************************************************************************
* Get From index
* overloaded function: allows individual cells to use the FREEZE stack where
* there is a time lag
* Arguments cell row and col
* Returns: i_from_stack when normal no freeze behaviour is going on
* i_freeze_stack when freezing is going on in this cell
* The GDM layers need to freeze to accumulate sufficient change to force a
* species shuffle or similar within a cell..this may happen with short time
* steps
****************************************************************************
int TRaft::GetFromIndex(const long int arg_i_row,
                        const long int arg_i_col)
{
  if(Map[arg_i_row][arg_i_col]->IsFrozen())
    return K_FREEZE;
  else
    return i_from_stack;
} // end overloaded Get From Index
/************************************************************************
* Returns pointer to the specified relative cell
* Arguments:
* const int arg_i_rel,   // number o fthe relative cell
* const int arg_i_row,   // cell which we are relative to
* const int arg_i_col    // cell which we are relative to
* Returns NULL if sea
************************************************************************
TGridCell_ptr TRaft::GetRelCell(const int arg_i_rel,
                                const long int arg_i_row,
                                const long int arg_i_col)
{
int i_r, i_c;

  i_r=s_rel_row[arg_i_rel]+ arg_i_row;
  if(i_r<0||i_r>=n_rows)
    return NULL; // OFF THE TOP OR BOTTOM OF THE GRID
                 // ASSUME OUT OF BOUNDS FOR TRANSVERSE MERCATOR

  i_c=s_rel_col[arg_i_rel]+ arg_i_col;
  // FOR COLUMNS WE HAVE TO WRAP
  if(i_c<0)
    {
    if(Run_Parameters->b_wrap)
      i_c+=n_cols;
    else
      return NULL;
    }  // end if wrapping off western edge
  else
    {
    if (i_c>=n_cols)
      {
      if(Run_Parameters->b_wrap)
        i_c-=n_cols;
      else
        return NULL;
      }// end if wrapping off Eastern edge
    }
return(Map[i_r][i_c]);
} // end Get Rel Cell

/************************************************************************
* KM Oct 2012
* Returns the row index of the specified relative cell
* Arguments:
* const int arg_i_rel,   // number of the relative cell
* const int arg_i_row,   // cell which we are relative to
* const int arg_i_col    // cell which we are relative to
* Returns NULL if sea
************************************************************************
int TRaft::GetRowRelCell(const int arg_i_rel,             // KM Oct 2012
                         const long int arg_i_row,
                         const long int arg_i_col)
{
int i_r, i_c;

  i_r=s_rel_row[arg_i_rel]+ arg_i_row;
  if(i_r<0||i_r>=n_rows)
    return NULL; // OFF THE TOP OR BOTTOM OF THE GRID
                 // ASSUME OUT OF BOUNDS FOR TRANSVERSE MERCATOR

  i_c=s_rel_col[arg_i_rel]+ arg_i_col;
  // FOR COLUMNS WE HAVE TO WRAP
  if(i_c<0)
    {
    if(Run_Parameters->b_wrap)
      i_c+=n_cols;
    else
      return NULL;
    }  // end if wrapping off western edge
  else
    {
    if (i_c>=n_cols)
      {
      if(Run_Parameters->b_wrap)
        i_c-=n_cols;
      else
        return NULL;
      }// end if wrapping off Eastern edge
    }
return(i_r);
} // end GetRowRelCell

/************************************************************************
* KM Oct 2012
* Returns the col index of the specified relative cell
* Arguments:
* const int arg_i_rel,   // number of the relative cell
* const int arg_i_row,   // cell which we are relative to
* const int arg_i_col    // cell which we are relative to
* Returns NULL if sea
************************************************************************
int TRaft::GetColRelCell(const int arg_i_rel,             // KM Oct 2012
                         const long int arg_i_row,
                         const long int arg_i_col)
{
int i_r, i_c;

  i_r=s_rel_row[arg_i_rel]+ arg_i_row;
  if(i_r<0||i_r>=n_rows)
    return NULL; // OFF THE TOP OR BOTTOM OF THE GRID
                 // ASSUME OUT OF BOUNDS FOR TRANSVERSE MERCATOR

  i_c=s_rel_col[arg_i_rel]+ arg_i_col;
  // FOR COLUMNS WE HAVE TO WRAP
  if(i_c<0)
    {
    if(Run_Parameters->b_wrap)
      i_c+=n_cols;
    else
      return NULL;
    }  // end if wrapping off western edge
  else
    {
    if (i_c>=n_cols)
      {
      if(Run_Parameters->b_wrap)
        i_c-=n_cols;
      else
        return NULL;
      }// end if wrapping off Eastern edge
    }
return(i_c);
} // end GetColRelCell

//*****************************************************************************
//* GET SPECIES LIST
//* Gets the actual species data from the specified cell
//* Arguments: row and col indices.
//* Note this can also be done FASTER by a direct call to Get Species Data
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
short* TRaft::GetSpeciesList(const long int arg_i_row,
                             const long int arg_i_col)
{
  return Map[arg_i_row][arg_i_col]->GetSpeciesData();
}
//*****************************************************************************
//* GET REL SPECIES LIST
//* Gets the actual species data from the specified relative cell
//* Arguments:
//* Relative index
//* row and col indices of the central cell, index o.
//* Note this can also be done FASTER by a direct call to Get Species Data
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
short* TRaft::GetRelSpeciesList( const int arg_i_rel,
                                 const long int arg_i_row,
                                 const long int arg_i_col)
{
   return GetRelCell(arg_i_rel,
                     arg_i_row,
                     arg_i_col)->GetSpeciesData();
}
//*****************************************************************************
//* Make Land
//* Uses the supplied file to allocate TGridCells to the pointers on the raft,
//* building a surface of land, but leaving the other pointers NULL
//* The file must be in the semi RLE encoded *.MCC format, with RLE for the NODATA
//* Any valid condition file will do, Pristine is good..but use the
//* From Condition grid by default
//* Arguments: string arg_s_condition_file  the full path to load from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::MakeLand(string arg_s_condition_file)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;   // stream to read from
int i;
long int i_col;
long int i_row;
long int i_land_count;  // to make species row reading easier
float f_buffer;         // reading stuff
float* f_ptr;
short* i_s_ptr;
short i_s;

  f_ptr=&(f_buffer);
  i_s_ptr=&(i_s);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  i_land_count=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
    {
	// Read run length and type
    InFileStream.read((char*)i_s_ptr,sizeof(short));
	if(i_s<0)
	  {
	  if(i_s<=K_EOF) // end file
        break;
      // nodata run
	  i_s=0-i_s; // make positive
      i_col+=i_s;

      if(i_col==n_cols-1)
        {
        // end of row
		n_land_in_row[i_row]=i_land_count;
        i_land_count=0; // reset land count for the next row
		i_row++;
        i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out layer by layer
      for(i=0;i<i_s;i++)
		{
        InFileStream.read((char*)f_ptr,sizeof(float)); //ignore this time
		i_col++;
		Map[i_row][i_col]=new TGridCell();
		// count valid cells
		i_land_count++;
        }  // end for i
      if(i_col==n_cols-1)
        {
        // end of row
        n_land_in_row[i_row]=i_land_count;
        i_land_count=0; // reset land count for the next row
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

  return b_success;
} // end func MAKE LAND

/*****************************************************************************
* CalculateRelativeRadialNeighbourhood
*  Calculates the relative position of the neighbours lying within a given
*  radius and stores as part of the grid
*  we only need to do this once, then we know the relative positions of all cells
*  within a neighbourhood
*  Argument: const float arg_f_radius: the radius IN grid units (KM) for the
****************************************************************************
void TRaft::CalculateRelativeRadialNeighbourhood(const float arg_f_radius)
{
int i_rad;     // radius
float f_rad;   // radius
int i_r, i_c;
int i_count;   // num cells count
int i_cells; // number of reauired relaitve cells
float f_dist;   // separation
double d_sq;
  if(Run_Parameters->i_units)      // metres.. radius is supplied in km
    f_rad=arg_f_radius*1000/Run_Parameters->f_cellsize ; // convert to number of cells
  else
    f_rad=arg_f_radius/Run_Parameters->f_cellsize ;

  i_rad=(int)f_rad;
  if(f_rad-i_rad>=0.5)
    i_rad++;
  // Calculate required size of array (with a little safety ring) and allocate
  i_cells=(int) PI*1.1*i_rad*i_rad;
  s_rel_row=AllocateShort1DArray(i_cells);
  s_rel_col=AllocateShort1DArray(i_cells);
  f_rel_dist=AllocateFloat1DArray(i_cells);
  f_rel_bearing=AllocateFloat1DArray(i_cells);
  d_rel_prob=AllocateDouble1DArray(i_cells);
  i_count=0;

  if(i_rad>0)
    {
    for(i_r=-i_rad; i_r<= i_rad; i_r++)
      {
      for(i_c=-i_rad; i_c<= i_rad; i_c++)
        {
        d_sq=double(i_r*i_r + i_c*i_c);
        if(d_sq>0)
          f_dist= sqrt(d_sq);
        else
          f_dist=0;
        if(f_dist <=f_rad)
          {
          // in radius add to list
          s_rel_row[i_count]=i_r;
          s_rel_col[i_count]=i_c;
          f_rel_dist[i_count]=f_dist;
          f_rel_bearing[i_count]=CalculateBearingToRel(i_r,
                                                       i_c);
          d_rel_prob[i_count]=CalculateDispersal(f_rel_dist[i_count]);
          i_count++; // first time through = 1
          }// end if in radius
        } // end for i_c
      } // end for i_r
    }// end if i_rad>0
  else
    {
    // always survey self
    s_rel_row[i_count]=0;
    s_rel_col[i_count]=0;
    i_count++; // first time through = 1
    }
  n_rel=i_count;
}

/*****************************************************************************
*  Load Relative Neighbourhood
*  Loads the relative position of the neighbours from a sample file
*  Argument: sample file *msf open, with TEXT format
*       i_num_samples/n
*       i_rel_row[0] i_rel_col[0] /n
*       i_rel_row[1] i_rel_col[1] /n    etc.
****************************************************************************
void TRaft::LoadRelativeNeighbourhood(string arg_s_sample_file)
{
ifstream InFile;
int i_count;
int n_samples; // reading only
long int i_r,i_c;   // buffers and easy locals
long int i_min_r, i_max_r;
double d_sq;
double d_prob_sum;

i_min_r=0;
i_max_r=0;

  // two passes so we can allocate nicely
  InFile.open(arg_s_sample_file.c_str());
  // Calculate required size of array (with a little safety ring) and allocate
  InFile>>n_samples;
  n_rel=n_samples;
  s_rel_row=AllocateShort1DArray(n_samples+1);
  s_rel_col=AllocateShort1DArray(n_samples+1);
  f_rel_dist=AllocateFloat1DArray(n_samples+1);
  f_rel_bearing=AllocateFloat1DArray(n_samples+1);
  d_rel_prob=AllocateDouble1DArray(n_samples+1);
  //now read in files. rel_row. space rel_col return
  i_count=0;
  d_prob_sum=0;
  for(i_count=0;i_count<n_samples;i_count++)
    {
    InFile>>i_r;
    s_rel_row[i_count]=i_r;
    InFile>>i_c;
    s_rel_col[i_count]=i_c;
    d_sq=double(i_r*i_r + i_c*i_c);
    if(d_sq>0)
      {
      f_rel_dist[i_count]= sqrt(d_sq);

      f_rel_bearing[i_count]=CalculateBearingToRel(i_r,
                                                   i_c);
      d_rel_prob[i_count]=CalculateDispersal(f_rel_dist[i_count]);
	  d_prob_sum+=d_rel_prob[i_count];
      }
    else
      {
      f_rel_dist[i_count]=0;
      f_rel_bearing[i_count]=0;
      d_rel_prob[i_count]=0;
      }
    if(i_r<i_min_r)
      i_min_r=i_r;
    if(i_r>i_max_r)
      i_max_r=i_r;
    }// end while i_count<n_samples
 for(i_count=0;i_count<n_samples;i_count++)
   d_rel_prob[i_count]/=d_prob_sum;                                             // is "/=" right?
 // we know the maximum range (e.g. -5 to + 5 = 10 (+ central 1 )= 11
 i_buffer_diameter=i_max_r-i_min_r+1;                                           //i_buffer_diameter=i_max_r-i_min_r+1; // store for future use
 i_buffer_radius=(i_buffer_diameter/2) + 1;                                     //i_buffer_radius= (i_max_r-i_min_r)/2;
}// end func load relative neighbourhood

/*****************************************************************************
* CalculateBearingToRel
*  Calculates the bearing in RADIANS to the relative position of the supplied
*  Relative cell
* A bit clunky with the ifs required so have taken it out of Calc Rel Neighbourhood
*  we only need to do this once, then we know the relative positions of all cells
*  within a neighbourhood
*  Argument: const float arg_f_radius: the radius IN grid units (KM) for the
* Returns FLOAT: the bearing in RADIANS
****************************************************************************
float TRaft::CalculateBearingToRel(const long int arg_i_rel_row,
                                        const long int arg_i_rel_col)
{
int i_r, i_c;
float f_bear;  // bearing in radians
  i_r=arg_i_rel_row;
  i_c= arg_i_rel_col;

  // Work out bearing
  // ROW MINUS NORTH, PLUS SOUTH
  // COL MINUS WEST, PLUS EAST
  //-row+col 0-0.5PI  N E
  //+row+col 0.5PI - 1PI  S E
  //+row-col 1PI-1.5 PI  S W
  //-row-col 1.5PI-2PI  N W
  if(i_r<0)
   {
   if(i_c<0)
     {
     //-row-col 1.5PI-2PI NORTH WEST
     f_bear=TWO_PI-atan(-(double)i_c/-(double)i_r);
     } // end if i_c<0)
   else
     {
     if (i_c>0)
       {
       //-row+col 0- 0.5PI  NORTH EAST
       f_bear=atan((double)i_c/-(double)i_r);
       } // end if i_c>0
     else
       {
       if(i_c==0)
         {
         f_bear=0;  //NORTH
         } // end if i_c==0
       } //end else
     } // end else ic

   } // end if i_r<0
  else
   {
   if(i_r>0)
     {
     if(i_c<0)
       {
       //+row-col 1PI-1.5 PI  SOUTH WEST
        f_bear=PI+atan(-(double)i_c/(double)i_r);
       } // end if i_c<0)
     else
       {
       if (i_c>0)
         {
         //+row+col 0.5PI - 1PI  SOUTH EAST
         f_bear=PI-atan((double)i_c/(double)i_r);
         } // end if i_c>0
       else
         {
         if(i_c==0)
           {
           f_bear=PI;  //SOUTH
           } // end if i_c==0
         } //end else
       } // end else ic
     } // end if i_r>0
   else
     {
     if(i_r==0)
       {
       if(i_c<0)
         f_bear=1.5*PI; // WEST
       else
         {
         if (i_c>0)
           f_bear=PI/2; // EAST
         else
           {
           if(i_c==0)
             f_bear=-1; //Self
           } //end else
         } // end else ic
       }// end if i_r==0
     }  // end else
   }// end else

  // Now simplify for now! 1 NNE 2 ENE 3 ESE 4 SSE 5 SSW 6 WSW 7 WNW 8 NNW
  // OR: 0 NNE 1 ENE 2 ESE 3 SSE 0 SSW 1 WSW 2 WNW 3 NNW

  if(f_bear<0.392699082)
    f_bear=0;       // N
  else
    {
    if(f_bear<1.178097245)
      f_bear=1;     //NE
    else
      {
      if(f_bear<1.963495408)
        f_bear=2;   //E
      else
        {
        if(f_bear<2.748893572)
          f_bear=3; //SE
        else
          {
          if(f_bear<3.534291735)
            f_bear=0;  //S
          else
            {
            if(f_bear<4.319689899)
              f_bear=1;   //SW
            else
              {
              if(f_bear<5.105088062)
                f_bear=2;   //W
              else
                {
                if(f_bear<5.890486225)
                  f_bear=3;   //NW
                else
                  {
                  f_bear=0; // N
                  }
                }
              }
            }
          }
        }
      }
    }
  return f_bear;
}
/*****************************************************************************
* CalculateDispersal
*  Calculates the probability of dispersal to a point arg_f_dist
* Returns FLOAT: the bearing in RADIANS
****************************************************************************
double TRaft::CalculateDispersal(const float arg_f_dist)
{
double d_prob;
double d_K, d_lambda,d_dist;

// Work entirely in metres
  if(Run_Parameters->i_units)
    {  // 1
    // metres
    // dist is "cellsize in metres* dist in cells"
    d_dist=Run_Parameters->f_cellsize*arg_f_dist;
    }
  else
    {  // 0
    // Kilometres
    // dist is "cellsize in kilometres* dist in cells*1000" to upgrade from m to km
    d_dist=Run_Parameters->f_cellsize*arg_f_dist*1000;
    }

  if(Run_Parameters->b_dispersal)
    {
    d_lambda=Run_Parameters->f_lambda;
    d_K=Run_Parameters->f_K;
    if(Run_Parameters->b_cauchy)
      {
      // Half cauchy
      if(d_dist>0)
        d_prob=2*d_K/((2*PI*d_dist)*PI*d_lambda*(1+(d_dist/d_lambda)*(d_dist/d_lambda)));
      else
        d_prob=1;
      }
    else
      {
      // Exponential
      if(d_dist>0)
        d_prob=d_K/( (2*PI*d_dist)*d_lambda*exp(d_dist/d_lambda));
      else
        d_prob=1;
      } //end if exp
    } // end if dispersal
  else
    d_prob=1;
  return d_prob;
}  // end fun Calc dispersal



//*****************************************************************************
//* Advance Time
//*
//* Moves the little sub time step
//* NO ALLOCATION or DESTRUCTION
//*
//*****************************************************************************
void TRaft::AdvanceTime(void)
{
long int i_row, i_col;
 i_step_count--;

 FlipStacks(); // make the current TO stack (e.g. 2020) into the FROM stack
 
 if(i_step_count==0)  // Reload next data level
   {
   RefreshAllSeries();
   }    // end if we need to load files
 else
   { // just add increments
   for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        Map[i_row][i_col]->IncrementGDMIncrement( i_to_stack,
                                                  n_layers);
        if(Run_Parameters->n_condition_steps>0)
          Map[i_row][i_col]->IncrementCondIncrement( i_to_stack);
        if(Run_Parameters->n_richness_steps>0)
          Map[i_row][i_col]->IncrementRichIncrement( i_to_stack);
        }
      } // end for i_col
    } // end for i_row
   } // end normal incrementing
}// end func Advance Time

/*****************************************************************************
* FlipStacks
*  makes i_from_stack=i_to_stack and vice versa
****************************************************************************
void TRaft::FlipStacks(void) // makes i_from_stack=i_to_stack dn vice versa
{
int i_stack;
  i_stack=i_from_stack;
  i_from_stack=i_to_stack;
  i_to_stack=i_stack;
}// end func FlipStacks

//*****************************************************************************
//* Initialise All Series
//* Loads the first two time points in the time series
//* Calculates the annual increment
//* Applies increment to current year to give present annual pair (THIS year NEXT year)
//* sets up the year counter
//* ( note this will be used by all the incrementing grids so we need to do all
//* the initialisings at the same time)
//* gets parameters from the RunParameters .
//* Doesn't do much if step is 0 (i.e. no increment stuff)
//*****************************************************************************
void TRaft::InitialiseAllSeries(void)
{
 i_from_stack=0;
 i_to_stack=1;
 i_step_years=Run_Parameters->n_sub_steps;
 i_step_count=i_step_years; //counts down for every advance

 // Starting point
 InitialiseGDMSeries();
 if(Run_Parameters->n_condition_steps>0)
   InitialiseCondSeries();
 if(Run_Parameters->n_richness_steps>0)
   InitialiseRichSeries();

 // Now the from stack and the to stack should be the first two years (e.g. 2010, 2011)
 // and i_step_count will be i_step_years  (e.g. 10)
 i_data_step=1;  //last step loaded
} // end func initialise All series
//*****************************************************************************
//* Refresh All Series
//* every i_step_years (e.g every 10 years)
//* Flips stacks, Loads all new data and recaculates the increment and  step counter
//*****************************************************************************
void TRaft::RefreshAllSeries(void)
{
string s_file;
 // Advance data step counter
 i_data_step++;  // index of the next FROM  (e.g. last one was 1 (from 0&1) now 1 (TO) becomes FROM and we load i_data_step(1) + 1 = 2
 //

 RefreshGDMSeries();
 if(Run_Parameters->n_condition_steps>0)
   RefreshCondSeries();
 if(Run_Parameters->n_richness_steps>0)
   RefreshRichSeries();

 i_step_count=i_step_years; // reset.. count down every year
} // end func refresh ALL series



//*****************************************************************************
//* Load Whole GDMStack
//* Loads the supplied file to its position on the raft, according to the
//* which index
//* The file must be in the semi RLE encoded *.MCG format, with RLE for the NODATA
//* ALLOCATES MEMORY in each TGridCell
//* Arguments:
//       const int arg_i_which is the current from or to stack
//       string arg_s_file  the full path to load the FLOAT MCG file from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::LoadWholeGDMStack(const int arg_i_which,  //K_FROM_STACK 0 and K_TO_STACK 1
                                  string arg_s_file)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
long int i_col;
long int i_row;
int i_layer;            // counter
float* f_ptr;   // reading pointer
float f_buffer; // place to read stuff we don't care about.

  i_s_ptr=&(i_s);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
    {
    // Read run length and type
    InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
      if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
      i_col+=i_s;

      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out layer by layer
      for(i=0;i<i_s;i++)
        {
        i_col++;
        if(Map[i_row][i_col]==NULL)
          {
          // NO DATA CELL HERE>>read out the file data to skip this cell
          for(i_layer=0;i_layer<n_layers;i_layer++)
            {
            f_ptr=&f_buffer;
            InFileStream.read((char*)f_ptr,sizeof(float));
            }
          }
        else
          {
           // Read in all the layers direct to the cell    ALLOCATE HERE
          Map[i_row][i_col]->LoadGDMStack(arg_i_which,     //Assume all is well and just set
                                          &InFileStream,
                                          n_layers);
          }
        }  // end for i
      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

return b_success;
}// end func Load whole GDMStack

//*****************************************************************************
//* Read In Whole GDMStack
//* Loads the supplied file to its position on the raft, according to the
//* which index
//* The file must be in the semi RLE encoded *.MCG format, with RLE for the NODATA
//*Deosn't allocate--just overwrites in each TGridCell
//* Arguments:
//       const int arg_i_which  current from and to stacks
//       string arg_s_file  the full path to load the FLOAT MCG file from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::ReadInWholeGDMStack(const int arg_i_which,  //K_FROM_STACK 0 and K_TO_STACK 1
                                  string arg_s_file)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
long int i_col;
long int i_row;
int i_layer;            // counter
float* f_ptr;   // reading pointer
float f_buffer; // place to read stuff we don't care about.

  i_s_ptr=&(i_s);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
    {
    // Read run length and type
    InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
      if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
      i_col+=i_s;

      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out layer by layer
      for(i=0;i<i_s;i++)
        {
        i_col++;
        if(Map[i_row][i_col]==NULL)
          {
          // NO DATA CELL HERE>>read out the file data to skip this cell
          for(i_layer=0;i_layer<n_layers;i_layer++)
            {
            f_ptr=&f_buffer;
            InFileStream.read((char*)f_ptr,sizeof(float));
            }
          }
        else
          {
          // Read in all the layers direct to the cell  NO ALLOC
		  Map[i_row][i_col]->SetGDMStack(arg_i_which,     //Assume all is well and just set
                                         &InFileStream,
                                         n_layers);
          }
        }  // end for i
      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

return b_success;
}// end func Read In whole GDMStack


//*****************************************************************************
//* Initialise GDM Series
//* Loads the first two time points in the GDM series
//* Calculates the annual increment
//* Applies increment to current year to give present annual pair (THIS year NEXT year)
//* sets up the year counter
//* ( note this will be used by all the incrementing grids so we need to do all
//* the initialisings at the same time)
//*****************************************************************************
// OLD VERSION //
/*
void TRaft::InitialiseGDMSeries(void)
{
string s_file;

 AllocateGDMStacks();
 MeltAllGDMStacks();
 // Starting point
 s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_gdm_stack[0]+".mcg";
 ReadInWholeGDMStack(i_from_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
                       s_file);
 //LoadWholeGDMStack(i_from_stack,      // Allocates here
 //                  s_file);
 // Next known data
 s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_gdm_stack[1]+".mcg";

 ReadInWholeGDMStack(i_to_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
                       s_file);
 //LoadWholeGDMStack(i_to_stack,
 //                  s_file);

 //Calculate the increment somewhere.. either on the raft or in the cell and save
 // create space to hold increment
 if(Run_Parameters->n_sub_steps>0)
   {
   AllocateGDMIncrement();
   CalculateGDMIncrement();
   }

} // end func initialise GDM series
*/
// NEW VERSION // NEW VERSION // 24 Oct 2012 // For single GDM stack //
/*
void TRaft::InitialiseGDMSeries(void)
{
string s_file;

 AllocateGDMStacks();
 // Starting point
 s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_gdm_stack[0]+".mcg";
 ReadInWholeGDMStack(i_from_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
					   s_file);

// if(Run_Parameters->n_sub_steps>0)




} // end func initialise GDM series
// END NEW VERSION //

//*****************************************************************************
//* Refresh GDM Series
//* every i_step_years (e.g every 10 years)
//* Flips stacks, Loads a new GDM stack, recaculates the increment and  step counter
//*****************************************************************************
void TRaft::RefreshGDMSeries(void)
{
string s_file;

 if(Run_Parameters->n_data_steps>1&&Run_Parameters->n_sub_steps>0)
   {
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_gdm_stack[i_data_step]+".mcg";
   ReadInWholeGDMStack(i_to_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
                       s_file);
   //we now need to recalculate the increment, and overwrite the TO stack
   if(Run_Parameters->n_sub_steps>0)
     CalculateGDMIncrement();   // This leaves the to stack in the correct state ( i.e. FROM stack + INCREMENT)
   }
} // end func refresh GDM series

/*****************************************************************************
* Allocate GDM Stacks.
* passes allocation call down to all grid cells
****************************************************************************/
// OLD VERSION //
/*
int TRaft::AllocateGDMStacks(void)
{
long int i_col;
long int i_row;
int i;    // stack
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        for(i=0;i<3;i++)                                                        //change to  "for(i=0;i<1;i++)"
          Map[i_row][i_col]->AllocateGDMStack(i,
                                              n_layers);
        }  // end if allocating
      } // end for i_col
    } // end for i_row
  return TRUE;
}// end allocate GDM stacks
*/
// NEW VERSION // NEW VERSION // 24 Oct 2012 // For single GDM stack //

/*int TRaft::AllocateGDMStacks(void)
{
long int i_col;
long int i_row;
int i;    // stack
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
		for(i=0;i<1;i++)                                                        //change to  "for(i=0;i<1;i++)"
		  Map[i_row][i_col]->AllocateGDMStack(i,
                                              n_layers);
        }  // end if allocating
	  } // end for i_col
    } // end for i_row
  return TRUE;
}// end allocate GDM stacks
// END NEW VERSION //
/*****************************************************************************
* Free GDM Stacks.
* passes free call down to all grid cells
****************************************************************************/
// OLD VERSION //
/*
void TRaft::FreeGDMStacks(void)
{
long int i_col;
long int i_row;
int i;    // stack
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        for(i=0;i<2;i++)                                                        //change to  "for(i=0;i<1;i++)"
          Map[i_row][i_col]->FreeGDMStack(i);
        }
      } // end for i_col
    } // end for i_row
}// end func free GDM Increment
*/
// NEW VERSION // NEW VERSION // 24 Oct 2012 // For single GDM stack //
/*
void TRaft::FreeGDMStacks(void)
{
long int i_col;
long int i_row;
int i;    // stack
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        for(i=0;i<1;i++)                                                        //change to  "for(i=0;i<1;i++)"
          Map[i_row][i_col]->FreeGDMStack(i);
        }
      } // end for i_col
    } // end for i_row
}// end func free GDM Increment

// END NEW VERSION //
/*****************************************************************************
* Mel All GDM Stacks.
* Sets all grid cells to K_MELT
****************************************************************************
void TRaft::MeltAllGDMStacks(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        Map[i_row][i_col]->SetFreeze(K_MELT,    // K_FREEZE or K_MELT
                                     0,  // Any old thing
                                     n_layers);
        }  // end if non null
      } // end for i_col
    } // end for i_row
}// end MeltAll GDM stacks
/*****************************************************************************
* Allocate GDM Increment.
* passes allocation call down to all grid cells
****************************************************************************
int TRaft::AllocateGDMIncrement(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        Map[i_row][i_col]->AllocateGDMIncrement(n_layers);
      } // end for i_col
    } // end for i_row
  return TRUE;
}// end allocate GDM increment
/*****************************************************************************
* Free GDM Increment.
* passes free call down to all grid cells
****************************************************************************
void TRaft::FreeGDMIncrement(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        Map[i_row][i_col]->FreeGDMIncrement();
      } // end for i_col
    } // end for i_row
}// end func free GDM Increment



//*****************************************************************************
//* Calculate GDM Increment
//*
//* 1. Calculates the annual increment between the From and To stacks, based on the
//* number of steps between
//* 2. Resets the To stack to be From+increment , so we're ready to roll
//* NO ALLOCATION or DESTRUCTION
//*
//*****************************************************************************
int TRaft::CalculateGDMIncrement(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        Map[i_row][i_col]->CalcGDMIncrement( Run_Parameters->n_sub_steps,
											 n_layers,
                                             i_from_stack);
        }
      } // end for i_col
    } // end for i_row
 // i_step_count=i_step_years; // reset.. count down every year
  return OK;
} // end func Calc GDM Increment


//*****************************************************************************
//* Load Condition Grid
//* Loads the supplied file to its position on the raft, according to the
//* which index
//* The file must be in the semi RLE encoded *.MCC format, with RLE for the NODATA
//* Arguments:
//       const int arg_i_which  use i_from_stack, and i_to_stack
//                       ZERO for the From Condition, ONE for the To Condition
//       string arg_s_condition_file  the full path to load from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::LoadConditionGrid(const int arg_i_which,
                                  string arg_s_condition_file)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
long int i_col;
long int i_row;
float f_cond;
float* f_ptr;

  i_s_ptr=&(i_s);
  f_ptr= &(f_cond);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_condition_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
    {
    // Read run length and type
    InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
      if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
      i_col+=i_s;

      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out
      for(i=0;i<i_s;i++)
        {
        InFileStream.read((char*)f_ptr,sizeof(float));
        i_col++;
        if(Map[i_row][i_col]!=NULL)      // only copy if there's land
          Map[i_row][i_col]->SetCondition(arg_i_which,
                                          f_cond);
        }  // end for i
      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

  return b_success;
} // end func Load Condition Grid

//*****************************************************************************
//* Initialise Cond Series
//* Loads the first two time points in the Cond series
//* Calculates the annual increment
//* Applies increment to current year to give present annual pair (THIS year NEXT year)
//* sets up the year counter
//* ( note this will be used by all the incrementing grids so we need to do all
//* the initialisings at the same time)
//*****************************************************************************
void TRaft::InitialiseCondSeries(void)
{
string s_file;

 // Starting point
 s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_condition[0]+".mcc";
 LoadConditionGrid(i_from_stack,
                   s_file);
 //Catch single condition option.. all the same
 if(Run_Parameters->n_condition_steps==1)
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_condition[0]+".mcc";
 else
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_condition[1]+".mcc";
 LoadConditionGrid(i_to_stack,
                   s_file);

 //Calculate the increment somewhere..
 if(Run_Parameters->n_sub_steps>0)
   CalculateCondIncrement();

} // end func InitialiseCondSeries

//*****************************************************************************
//* Refresh Cond Series
//* every i_step_years (e.g every 10 years)
//* Flips stacks, Loads a new GDM stack, recaculates the increment and  step counter
//*****************************************************************************
void TRaft::RefreshCondSeries(void)
{
string s_file;
 if(Run_Parameters->n_condition_steps>1&&Run_Parameters->n_sub_steps>0)
   {
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_condition[i_data_step]+".mcc";
   LoadConditionGrid(i_to_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
                     s_file);
   //we now need to recalculate the increment, and overwrite the TO stack
   if(Run_Parameters->n_sub_steps>0)
     CalculateCondIncrement();   // This leaves the to stack in the correct state ( i.e. FROM stack + INCREMENT)
   }
} // end func Refresh Cond series

//*****************************************************************************
//* Calculate Condition Increment
//*
//* 1. Calculates the annual increment between the From and To stacks, based on the
//* number of steps between
//*
//* NO ALLOCATION or DESTRUCTION
//*
//*****************************************************************************
int TRaft::CalculateCondIncrement(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        Map[i_row][i_col]->CalcCondIncrement( Run_Parameters->n_sub_steps,
                                              i_from_stack);
        }
      } // end for i_col
    } // end for i_row
  return OK;
} // end func Calc Cond Increment
//*****************************************************************************
//* Load Richness Grid
//* Loads the supplied file to its position on the raft, according to the
//* which index
//* The file must be in the semi RLE encoded *.MCC format, with RLE for the NODATA
//* Arguments:
//       const int arg_i_which  use i_from_stack, and i_to_stack
//                       ZERO for the From Condition, ONE for the To Condition
//       string arg_s_richness_file  the full path to load from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::LoadRichnessGrid(const int arg_i_which,
                                  string arg_s_richness_file)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
float f_rich;
float* f_ptr;
long int i_col;
long int i_row;

  i_s_ptr=&(i_s);
  f_ptr=&f_rich;
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_richness_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
    {
    // Read run length and type
    InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
      if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
      i_col+=i_s;

      if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out
      for(i=0;i<i_s;i++)
        {
        InFileStream.read((char*)f_ptr,sizeof(float));
        i_col++;
        //TEMP - don't allow zero richness
        if(f_rich<1)
          f_rich = 1;
        //END TEMP
        if(Map[i_row][i_col]!=NULL) // only copy to Raft if this cell is land
		  Map[i_row][i_col]->SetRichness(arg_i_which,
                                         f_rich);
        }  // end for i
      if(i_col==n_cols-1)
		{
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

  return b_success;
} // end func Load Richness Grid

//*****************************************************************************
//* Initialise Rich Series
//* Loads the first two time points in the Rich series
//* Calculates the annual increment
//* Applies increment to current year to give present annual pair (THIS year NEXT year)
//* sets up the year counter
//* ( note this will be used by all the incrementing grids so we need to do all
//* the initialisings at the same time)
//*****************************************************************************
void TRaft::InitialiseRichSeries(void)
{
string s_file;

 // Starting point
 s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_richness[0]+".mcr";
 LoadRichnessGrid(i_from_stack,
                   s_file);
 if(Run_Parameters->n_richness_steps==1)
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_richness[0]+".mcr";
 else
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_richness[1]+".mcr";
 LoadRichnessGrid(i_to_stack,
                   s_file);

 //Calculate the increment somewhere..
 if(Run_Parameters->n_sub_steps>0)
   CalculateRichIncrement();

} // end func InitialiseRichSeries

//*****************************************************************************
//* Refresh Rich Series
//* every i_step_years (e.g every 10 years)
//* Flips stacks, Loads a new serires, recaculates the increment and  step counter
//*****************************************************************************
void TRaft::RefreshRichSeries(void)
{
string s_file;

 if(Run_Parameters->n_richness_steps>1&&Run_Parameters->n_sub_steps>0)
   {
   s_file= Run_Parameters->s_input_dir+K_FOLDER+ Run_Parameters->s_richness[i_data_step]+".mcr";
   LoadRichnessGrid(i_to_stack,   // just overwrites current To stack 9 which is now holding the  From data from the previous time step
                     s_file);
   //we now need to recalculate the increment, and overwrite the TO stack
   if(Run_Parameters->n_sub_steps>0)
     CalculateRichIncrement();   // This leaves the to stack in the correct state ( i.e. FROM stack + INCREMENT)
   }
} // end func Refresh Rich series


//*****************************************************************************
//* Calculate Richness Increment
//*
//* 1. Calculates the annual increment between the From and To stacks, based on the
//* number of steps between
//*
//* NO ALLOCATION or DESTRUCTION
//*
//*****************************************************************************
int TRaft::CalculateRichIncrement(void)
{
long int i_col;
long int i_row;
  for(i_row=0;i_row<n_rows;i_row++)
    {
    for(i_col=0;i_col<n_cols;i_col++)
      {
      if(Map[i_row][i_col]!=NULL)
        {
        Map[i_row][i_col]->CalcRichIncrement( Run_Parameters->n_sub_steps,
                                              i_from_stack);
        }
      } // end for i_col
    } // end for i_row
  return OK;
} // end func Calc Rich Increment
//*****************************************************************************
//* InitialiseRealisedCompositionBuffer
//*
//* Stores file and sets up the species buffer to hold the first lot of data
//* if the start row is low, then the buffer will start at 0 and may be shallow
//*
//* ALLOCATES SPACE FOR SPECIES DATA
//* Argument; arg_i_start_row... the row at the centre of the buffer
//*
//*****************************************************************************
void TRaft::InitialiseRealisedCompositionBuffer(ifstream* arg_composition_file,
                                                 const long int arg_i_start_row)
{
long int i_row;
  // Pass composition file name on
  SetCompositionFileStream(arg_composition_file);
  // check our bounds.. probably i_min will be on 0
  i_min_row= arg_i_start_row- i_buffer_radius;
  if(i_min_row<0) i_min_row=0;
  i_max_row= arg_i_start_row+ i_buffer_radius;
  if(i_max_row>n_rows) i_max_row=n_rows;    // note this is like n_rows.. one more than the actual index

  // Load in the buffer
  for(i_row=i_min_row;i_row<i_max_row;i_row++)
    {
    LoadSpeciesRow(i_row);
    }
} // end func InitialiseRealisedCompositionBuffer
//*****************************************************************************
//* DestroyRealisedCompositionBuffer
//*
//* Frees up the rest of the buffer
//*
//* FREES SPACE FOR SPECIES DATA
//*
//*****************************************************************************
void TRaft::DestroyRealisedCompositionBuffer(void)
{
long int i_row;
  for(i_row=i_min_row;i_row<i_max_row;i_row++)
    {
    FreeSpeciesRow(i_row);
    }
} // end func DestroyRealisedCompositionBuffer
//*****************************************************************************
//* MoveRealisedCompositionBufferToNextRow
//*
//* Reads in the requested line of data and rolls on down
//*
//*
//*****************************************************************************
void TRaft::MoveRealisedCompositionBufferToNextRow(void)
{
long int i_new_max;
  //OK, we have to do two things..
  // 1. delete the min row and update min row to be one more
  // But only of we have reached the full buffer size
  if((i_max_row-i_min_row)>i_buffer_diameter)                                  // if((i_max_row-i_min_row)>=i_buffer_diameter)
    {
    FreeSpeciesRow(i_min_row);
    i_min_row++;
    }
  // 2. load a new row on the end.. if not possible just freeze
  i_new_max=i_max_row+1;
  if(i_new_max>n_rows)
    i_max_row=n_rows;  // run off edge of grid
  else
    {
    LoadSpeciesRow(i_max_row);    // note i_max_row is 1 more than the last row in the buffer
    i_max_row=i_new_max;
    }
} // end func MoveRealisedCompositionBufferToNextRow
//*****************************************************************************
//* Load Species Row
//*
//* Reads in the requested line of data from the current composition file
//* SETS the Richness for the current cell to the read in value
//* ALLOCATES SPACE FOR SPECIES DATA
//* Argument; arg_i_row... the row to load
//*
//*****************************************************************************
void TRaft::LoadSpeciesRow(const long int arg_i_row)
{// Declare local variables
short int* si_ptr;
short int si_buffer;    // short buffer to read into
long int i_row_index,i_col_index;
long int i_land_count;  // counter to simplify row reading
short int* si_spp_data;
int i_sp;
short int si_Richness;    // species richness
short int si_species[500]; // NOTE WE ARE CURRENTLY LIMITED TO 500 Species

si_ptr = &si_buffer;
  i_land_count=0;
  while(i_land_count<n_land_in_row[arg_i_row])
    {
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    i_row_index = si_buffer;
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    i_col_index = si_buffer;
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    si_Richness = si_buffer;

    for(i_sp=0; i_sp<si_Richness; i_sp++)
      {
      CompositionFileStream->read((char*)si_ptr,sizeof(short));
      si_species[i_sp] = si_buffer;    // typecast to short...not enough species to cause problems
      }  // end for present species
      
    // update the Richness for this cell to reflect the REALISED richness
    Map[i_row_index][i_col_index]->SetRichness(i_from_stack,     //Assume From or To
                                              (float)si_Richness);
    // At this point we know the coordinates of the point and
    // how many (s_richness) & which (s_species) species are at each point
    // create a suitable array  to hold the data
    si_spp_data= AllocateShort1DArray(si_Richness);
	Map[i_row_index][i_col_index]->SetSpeciesData(si_spp_data);
    // copy across to allocated array
    for(i_sp=0;i_sp<si_Richness;i_sp++)
      si_spp_data[i_sp]=si_species[i_sp];
    i_land_count++;
    } // end while i_land_count
} // end Load Sp Row
//*****************************************************************************
//* Load Known Species Data Row    // KM Oct 2012
//*
//* Reads in the requested line of data from the current composition file
//* SETS the Richness for the current cell to the read in value
//* ALLOCATES SPACE FOR SPECIES DATA - but only where there are known data present
//* Argument; arg_i_row... the row to load
//*
//*****************************************************************************
void TRaft::LoadKnownSpeciesDataRow(const long int arg_i_row)
{// Declare local variables
short int* si_ptr;
short int si_buffer;    // short buffer to read into
long int i_row_index,i_col_index;
long int i_land_count;  // counter to simplify row reading
short int* si_known_spp_data;
int i_sp;
short int si_Richness;    // species richness
short int si_species[1000]; // NOTE WE ARE CURRENTLY LIMITED TO 1000 Species

si_ptr = &si_buffer;
  i_land_count=0;
  while(i_land_count<n_land_in_row[arg_i_row])
    {
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    i_row_index = si_buffer;
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    i_col_index = si_buffer;
    CompositionFileStream->read((char*)si_ptr,sizeof(short));
    si_Richness = si_buffer;

    for(i_sp=0; i_sp<si_Richness; i_sp++)
      {
      CompositionFileStream->read((char*)si_ptr,sizeof(short));
      si_species[i_sp] = si_buffer;    // typecast to short...not enough species to cause problems
      }  // end for present species

    // update the Richness for this cell to reflect the REALISED richness
    Map[i_row_index][i_col_index]->SetRichness(i_to_stack,     //Assume To for known data
                                              (float)si_Richness);
    // At this point we know the coordinates of the point and
    // how many (s_richness) & which (s_species) species are at each point
    // create a suitable array  to hold the data
    if(si_Richness > 0)
      {
      si_known_spp_data= AllocateShort1DArray(si_Richness);
      Map[i_row_index][i_col_index]->SetKnownSpeciesData(si_known_spp_data);
      // copy across to allocated array
      for(i_sp=0;i_sp<si_Richness;i_sp++)
        si_known_spp_data[i_sp]=si_species[i_sp];
      } // end if si_Richness > 0

    i_land_count++;
    } // end while i_land_count
} // end LoadKnownSpeciesDataRow

//*****************************************************************************
//* Load Empty Species Data Row    // KM Oct 2012
//*
//* ALLOCATES SPACE FOR SPECIES DATA THAT ARE YET TO BE GENERATED, but where the
//* species richness is known
//* Argument; arg_i_row... the row to load
//*
//*****************************************************************************
void TRaft::LoadEmptySpeciesDataRow(const long int arg_i_row)
{// Declare local variables
long int i_row_index,i_col_index;
short int* si_spp_data;
int i_sp;
short int si_Richness;    // species richness
long int i_col;

  // Loop through all the columns in the row. If the species richness is > 0,
  // allocate a short array of the right length to catch the species data
  for(i_col=0;i_col<n_cols;i_col++)
	{
	if(Map[arg_i_row][i_col]!=NULL)
	  {
	  si_Richness = Map[arg_i_row][i_col]->GetRichness(K_FROM_STACK);
	  if(si_Richness > 0)
		{
		si_spp_data= AllocateShort1DArray(si_Richness);
		Map[arg_i_row][i_col]->SetSpeciesData(si_spp_data);
		// Fill in the allocated array with zeros
		for(i_sp=0;i_sp<si_Richness;i_sp++)
		  si_spp_data[i_sp]=0;
		} // end if si_Richness > 0
	  } // end if Map[i_row][i_col]!=NULL
	} // end for i_col

} // end LoadEmptySpeciesDataRow

//*****************************************************************************
//* Free Species Row
//*
//* Deletes the supplied row line of data from the Raft
//* Only deletes valid data in valid cells
//* ALLOCATES SPACE FOR SPECIES DATA
//* Argument; arg_i_row... the row to delete
//*
//*****************************************************************************
void TRaft::FreeSpeciesRow(const long int arg_i_row)
{
long int i_col;
short int* si_spp_data;

  for(i_col=0;i_col<n_cols;i_col++)
    {
    if(Map[arg_i_row][i_col]!=NULL)
      {
      si_spp_data= Map[arg_i_row][i_col]->GetSpeciesData();
      if(si_spp_data!=NULL)
        FreeShort1DArray(si_spp_data);
      si_spp_data=NULL;
      }
    }
} // end FreeSpeciesRow

//*****************************************************************************
//* Free Known Data Species Row    // KM Oct 2012
//*
//* Deletes the supplied row line of data from the Raft
//* Only deletes valid data in valid cells
//* ALLOCATES SPACE FOR SPECIES DATA
//* Argument; arg_i_row... the row to delete
//*
//*****************************************************************************
void TRaft::FreeKnownSpeciesDataRow(const long int arg_i_row)
{
long int i_col;
short int* si_known_spp_data;

  for(i_col=0;i_col<n_cols;i_col++)
    {
	if(Map[arg_i_row][i_col]!=NULL)
	  {
	  si_known_spp_data= Map[arg_i_row][i_col]->GetKnownSpeciesData();
	  if(si_known_spp_data!=NULL)
		FreeShort1DArray(si_known_spp_data);
	  si_known_spp_data=NULL;
	  }
    }
}  // end FreeKnownSpeciesDataRow

//---------------------------------------------------------------------------
//*****************************************************************************
//* Load Empty Raft Data    // KM May 2014
//*
//* ALLOCATES SPACE FOR DATA REQUIRED TO BE HELD ON THE RAFT
//* Argument; arg_i_row... the row to load
//*
//*****************************************************************************/
void TRaft::LoadEmptyRaftData(int arg_n_rows,
							  int arg_n_cols,
							  int arg_n_env_variables)
{
  // Declare local variables
  long int i_row, i_col, i_env;
  float* f_env_data_array;
  double* d_env_current_z_array;
  double* d_env_current_Vp_array;
  double* d_env_next_z_array;
  double* d_env_next_Vp_array;

  // Loop through the Raft, allocate the necessary grids where required, fill
  // them with zeros
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(Map[i_row][i_col]!=NULL)
		{
		f_env_data_array=AllocateFloat1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetEnvData(f_env_data_array);
		d_env_current_z_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetCurrentZ(d_env_current_z_array);
		d_env_current_Vp_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetCurrentVp(d_env_current_Vp_array);
		d_env_next_z_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetNextZ(d_env_next_z_array);
		d_env_next_Vp_array=AllocateDouble1DArray(arg_n_env_variables);
		Map[i_row][i_col]->SetNextVp(d_env_next_Vp_array);
		// Fill in the allocated arrays with zeros
		for(i_env=0;i_env<arg_n_env_variables;i_env++)
		  {
		  f_env_data_array[i_env]=0;
		  d_env_current_z_array[i_env]=0;
		  d_env_current_Vp_array[i_env]=0;
		  d_env_next_z_array[i_env]=0;
		  d_env_next_Vp_array[i_env]=0;
		  }
		} // end if Map[i_row][i_col]!=NULL
	  } // end for i_col
	} // end for i_row

} // end LoadEmptyRaftData
//*****************************************************************************
//*****************************************************************************
//* FREE Raft Data    // KM May 2014
//*
//* FREES DATA ALLOCATED TO THE RAFT
//* Argument;
//*
//*****************************************************************************
void TRaft::FreeRaftData(int arg_n_rows,
						 int arg_n_cols)
{
  // Declare local variables
  long int i_row, i_col;
  float* f_env_data_array;
  double* d_env_current_z_array;
  double* d_env_current_Vp_array;
  double* d_env_next_z_array;
  double* d_env_next_Vp_array;

  // Loop through the Raft, allocate the necessary grids where required, fill
  // them with zeros
  for(i_row=0;i_row<arg_n_rows;i_row++)
	{
	for(i_col=0;i_col<arg_n_cols;i_col++)
	  {
	  if(Map[i_row][i_col]!=NULL)
		{
		// The Env data array
		f_env_data_array= Map[i_row][i_col]->GetEnvData();
		if(f_env_data_array!=NULL)
		  FreeFloat1DArray(f_env_data_array);
		f_env_data_array=NULL;
		// The Env current z array
		d_env_current_z_array= Map[i_row][i_col]->GetCurrentZ();
		if(d_env_current_z_array!=NULL)
		  FreeDouble1DArray(d_env_current_z_array);
		d_env_current_z_array=NULL;
		// The Env current Vp array
		d_env_current_Vp_array= Map[i_row][i_col]->GetCurrentVp();
		if(d_env_current_Vp_array!=NULL)
		  FreeDouble1DArray(d_env_current_Vp_array);
		d_env_current_Vp_array=NULL;
		// The Env next z array
		d_env_next_z_array= Map[i_row][i_col]->GetNextZ();
		if(d_env_next_z_array!=NULL)
		  FreeDouble1DArray(d_env_next_z_array);
		d_env_next_z_array=NULL;
		// The Env next Vp array
		d_env_next_Vp_array= Map[i_row][i_col]->GetNextVp();
		if(d_env_next_Vp_array!=NULL)
		  FreeDouble1DArray(d_env_next_Vp_array);
		d_env_next_Vp_array=NULL;
		} // end if Map[i_row][i_col]!=NULL
	  } // end for i_col
	} // end for i_row

} // end FreeRaftData
//*****************************************************************************

//*****************************************************************************
//* Read In Whole ENVStack
//* Loads the supplied file to its position on the raft
//* The file must be in the semi RLE encoded *.MCG format, with RLE for the NODATA
//* Doesn't allocate--just overwrites in each TGridCell
//* Arguments:
//       string arg_s_file  the full path to load the FLOAT MCG file from
//* Returns: ONE success, ZERO Failure
//*****************************************************************************
int TRaft::ReadInWholeENVStack(string arg_s_file,
							   int arg_n_env_variables)
{
int b_success=FALSE;  // return flag
ifstream InFileStream;
int i;
short* i_s_ptr;
short i_s;
long int i_col;
long int i_row;
int i_layer;            // counter
float* f_ptr;   // reading pointer
float f_buffer; // place to read stuff we don't care about.

  i_s_ptr=&(i_s);
  i_col=MINUS_ONE;  // first increment =0
  i_row=0;
  // Open binary file up for binary reading
  InFileStream.open(arg_s_file.c_str(), ios::in |ios::binary);
  while(!InFileStream.eof())
	{
    // Read run length and type
	InFileStream.read((char*)i_s_ptr,sizeof(short));
    if(i_s<0)
      {
	  if(i_s<=K_EOF) // end file
        break;
      // nodata run
      i_s=0-i_s; // make positive
	  i_col+=i_s;

	  if(i_col==n_cols-1)
        {
		i_row++;
		i_col=MINUS_ONE;  // first increment =0
        }
      }
    else
      {
      // A run of data, read out layer by layer
	  for(i=0;i<i_s;i++)
		{
        i_col++;
		if(Map[i_row][i_col]==NULL)
          {
          // NO DATA CELL HERE>>read out the file data to skip this cell
		  for(i_layer=0;i_layer<arg_n_env_variables;i_layer++)
			{
			f_ptr=&f_buffer;
			InFileStream.read((char*)f_ptr,sizeof(float));
            }
          }
        else
		  {
          // Read in all the layers direct to the cell  NO ALLOC
		  Map[i_row][i_col]->SetENVStack(&InFileStream,
										 arg_n_env_variables);
          }
		}  // end for i
	  if(i_col==n_cols-1)
        {
        i_row++;
        i_col=MINUS_ONE; // first increment =0
        }
      } // end elsedata
    if(i_row==n_rows)
      break;
    } // end while eof
  //Close input file
  if(InFileStream.is_open())
    InFileStream.close();

return b_success;
}// end func Read In whole GDMStack


//*****************************************************************************

