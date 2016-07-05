//---------------------------------------------------------------------------

#ifndef CompactorFunctionsH
#define CompactorFunctionsH
using std::string;

// Float versions
int ConvertASCIIGridStackToCompactFloatBinaryFile(int i_nothing,
														string arg_str_parameter_filename);
/*
extern int ConvertConditionGridToCompactFloatBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);

extern int ConvertASCIIConditionGridToCompactFloatBinaryFile(string arg_s_condition_file,
                                                             string arg_s_res_path);
*/
int ConvertASCIIRichnessGridToCompactFloatBinaryFile(int i_nothing,
													  string arg_str_parameter_filename);

// Unsigned short versions
/*
extern int ConvertGridStackToCompactShortBinaryFile(const long int arg_i_num_rows,
                                                    const long int arg_i_num_cols,
                                                    const int arg_i_num_layers,
                                                    string arg_s_layer_files[],
                                                    string arg_s_res_path);

extern int ConvertConditionGridToCompactUShortBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);
extern int ConvertRichnessGridToCompactUShortBinaryFile(const long int arg_i_num_rows,
                                                        const long int arg_i_num_cols,
                                                        string arg_s_condition_file,
                                                        string arg_s_res_path);
*/
//---------------------------------------------------------------------------
#endif
