// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.

//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <Sift.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

struct kdtree
{
    int val, index;
    kdtree *l, *r;
    vector < SiftDescriptor > descriptors;
};

struct ans
{
    SiftDescriptor in_match;
    double val;
    double mined, min2ed;
    bool ismatch;
};

struct same
{
    int in_index, out_index;
    double val;
};

struct list
{
    string fname;
    int ct_match, ct_descriptors;
    float avg;
};

class compare
{
    public:
    	compare (int ind): i(ind) {}
    	bool operator()(const SiftDescriptor & a, const SiftDescriptor & b) {return (a.descriptor[i] < b.descriptor[i]);}
    private:
    	int i;
};

void gen_kdtree(vector < SiftDescriptor >, kdtree **, int, int);
void print_kdtree(kdtree **);

double euclidean_dist(SiftDescriptor &, SiftDescriptor &);
void search_kdtree(kdtree **, SiftDescriptor, double &, double, SiftDescriptor &);
int check_out_image(kdtree **, vector < SiftDescriptor >, double, vector < ans > &);

void feature_mark(CImg<double> &, vector < SiftDescriptor > &);
void feature_mark_append(CImg<double> , CImg<double> &, CImg<double> &, vector <ans> &, vector < SiftDescriptor > &);

CImg<double> matrixmult( CImg<double> , CImg<double> );
void warping(CImg<double> ,string );

vector < vector <float> > k_descriptorgenerator(vector < SiftDescriptor > &descriptors,int k,int w, vector < vector <double> > x)
{
    vector < vector <float> > k_descriptor;
    vector < float > temp_k_descriptor;
    for(int u=0; u<descriptors.size(); u++)
    {
        temp_k_descriptor.clear();
        for(int v=0; v<k; v++)
        {
            double xval = 0;
            for(int z=0; z<128; z++)
            {
                xval += (descriptors[u].descriptor[z] * x[v][z]);
            }
            xval /= w;
            xval = floor(xval);
            temp_k_descriptor.push_back(xval);
        }
        k_descriptor.push_back(temp_k_descriptor);
    }
    return k_descriptor;
}

vector < vector <double> > generaterandomvector(int k)
{
    CImg<double> r(k, 128);
    r.rand(0,1.0);

    vector < vector <double> > x;
    vector<double> tempx;
    for (int u=0; u<k; u++)
    {
        tempx.clear();
        for (int v=0; v<128; v++)
            tempx.push_back(r(v, u));
        x.push_back(tempx);
    }
    return x;
};

int main(int argc, char **argv)
{
    try 
    {
        if(argc < 3)
        {
            cout << "Insufficent number of arguments; correct usage:" << endl;
            cout << "    a2 part_id queryimage..." << endl;
            return -1;
        }
        
        string part = argv[1];
        string inputFile = argv[2];
    
        if(part == "part1" || part == "part2")
        {

            CImg<double> input_image(inputFile.c_str());
            CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
            vector < SiftDescriptor > descriptors = Sift::compute_sift(gray);
          
            vector < CImg<double> > output_images, output_gray;
            vector < vector < SiftDescriptor > > output_descriptors;
            for ( int of=3; of<argc; of++)
            {
                CImg<double> temp(argv[of]);
                output_images.push_back(temp);
                output_gray.push_back(output_images[of-3].get_RGBtoHSI().get_channel(2));
                output_descriptors.push_back(Sift::compute_sift(output_gray[of-3]));
            }
            
            vector < ans > answer;
            vector < CImg<double> > appended_images( output_images.size(), input_image);
            vector < vector <float> > in_k_descriptor, out_k_descriptor;

            string out_name, out_path = "./output-a2-images/part1_images/";
            if(part == "part2")
                out_path = "./output-a2-images/part2_images/";
            inputFile = inputFile.substr( inputFile.find_last_of("/")+1, string::npos);
            inputFile = inputFile.substr( 0, inputFile.find('.'));
            cout << "input file: " << inputFile << endl;

            int k=10, w=250;
            vector < vector <double> > x;
            if(part == "part2")
            {
                x = generaterandomvector(k);
                in_k_descriptor = k_descriptorgenerator(descriptors,k,w,x);
            }

            kdtree *root= new kdtree;
            int level = 10;
            gen_kdtree( descriptors, &root, 0, level);
                //print_kdtree( &root);
            
            int out_ct = output_descriptors.size(), count_match, j, sz, q[4];
            double thresh = 0.75, ratio;
            vector < list > desclist;
            list temp;

            for(int i=0; i<out_ct; i++)
            {
                count_match = 0;
                answer.clear();
                sz = output_descriptors[i].size();
                //cout << sz << " " << answer.size() << " " << count_match << endl;
                if(part == "part1")
                {
                    count_match = check_out_image( &root, output_descriptors[i], thresh, answer);
                  
                    feature_mark_append(input_image, output_images[i], appended_images[i], answer, output_descriptors[i]);
                }
                else
                {
                  
                out_k_descriptor = k_descriptorgenerator(output_descriptors[i],k,w,x);
                  
                answer.clear();
                count_match = 0;
                for(int a=0; a<out_k_descriptor.size(); a++)
                {
                    int b; 
                    for(b=0; b<in_k_descriptor.size(); b++)
                    {
                        int c=0;
                        for(; c<k; c++)
                        {
                            if(out_k_descriptor[a][c] != in_k_descriptor[b][c])
                                break;
                        }
                        if(c==k)
                          break;
                    }
                    ans temp;
                    if(b==in_k_descriptor.size())
                    {
                        temp.ismatch = false;
                        temp.val = -1;
                    }
                    else
                    {
                        double min, min2;
                        min2 = -1;
                        min = FLT_MAX;
                        SiftDescriptor close, close2;
                        search_kdtree( &root, output_descriptors[i][a], min, min2, close);
                        
                        min2 = FLT_MAX;
                        search_kdtree( &root, output_descriptors[i][a], min2, min, close2);
                        temp.mined = min;
                        temp.min2ed = min2;
                        temp.val = -1;
                        temp.ismatch = false;
                        min = min/min2;
                        if( (close.row == descriptors[b].row && close.col == descriptors[b].col) || (close2.row == descriptors[b].row && close2.col == descriptors[b].col) )
                        {
                            //if(min<=thresh)
                            //{    
                                temp.ismatch = true;
                                temp.in_match = descriptors[b];
                                temp.val = min;
                                count_match++;
                            //}
                        }
                        int d;
                        for(d=0; d<answer.size() && temp.ismatch; d++)
                            if(answer[d].in_match.row == temp.in_match.row && answer[d].in_match.col == temp.in_match.col && answer[d].ismatch)
                                break;
                        if(d < answer.size() && temp.ismatch && answer[d].ismatch)
                        {
                            if(answer[d].mined > temp.mined)
                            {    
                              answer[d].ismatch = false;
                              count_match--;
                            }
                            else
                            {  
                              temp.ismatch = false;
                              count_match--;
                            }
                        }
                    }
                    answer.push_back(temp);
                }
                //cout << "count_match: " << count_match << endl;
                
                CImg<double> A(8,8), B(1,8), h(3,3), newin(1,3), best_h(3,3);
		            best_h(0,0)=1;
		            best_h(1,0)=0;
		            best_h(2,0)=0;
		            best_h(0,1)=0;
		            best_h(1,1)=1;
		            best_h(2,1)=0;
		            best_h(0,2)=0;
		            best_h(1,2)=0;
		            best_h(2,2)=1;
		            h(2,2)=1;
		            newin(0,2)=1;
		            double dist, no, best_no=FLT_MAX;
		            for(int u=0; u<400 && count_match >=4 ; u++)
		            {
		                do
		                {
		                    q[0] = rand() % sz;
		                } while (!answer[q[0]].ismatch);
		                do
		                {
		                    q[1] = rand() % sz;
		                }while(q[1] == q[0] || !answer[q[1]].ismatch);
		                do
		                {
		                    q[2] = rand() % sz;
		                }while(q[2] == q[0] || q[2] == q[1] || !answer[q[2]].ismatch);
		                do
		                {
		                    q[3] = rand() % sz;
		                }while(q[3] == q[0] || q[3] == q[1] || q[3] == q[2] || !answer[q[3]].ismatch);
		                for( int r=0; r<8; r++)
		                {
		                    if(r%2 == 0)
		                    {
		                        A(0,r)=answer[q[r/2]].in_match.col;
		                        A(1,r)=answer[q[r/2]].in_match.row;
		                        A(2,r)=1;
		                        for(int t=3; t<6; ++t)
		                            A(t,r)=0;
		                        double x = -1 * output_descriptors[i][q[r/2]].col;
		                        A(6,r)=(x * answer[q[r/2]].in_match.col);
		                        A(7,r)=(x * answer[q[r/2]].in_match.row);
		                        B(0, r) = output_descriptors[i][q[r/2]].col;
		                        B(0, r+1) = output_descriptors[i][q[r/2]].row;
		                    }
		                    else
		                    {   
		                        for(int t=0; t<3; ++t)
		                            A(t,r)=0;
		                        A(3,r)=answer[q[r/2]].in_match.col;
		                        A(4,r)=answer[q[r/2]].in_match.row;
		                        A(5,r)=1;
		                        double y = -1 * output_descriptors[i][q[r/2]].row;
		                        A(6,r)=(y * answer[q[r/2]].in_match.col);
		                        A(7,r)=(y * answer[q[r/2]].in_match.row);
		                    }
		                }
		                CImg<double> tA = A.invert();
		                CImg<double> C=tA*B;
		                for(int tr=0; tr<8; tr++)
		                  h(tr%3,tr/3)=C(0,tr);
		                no=0;
		                for(int e=0; e<sz; e++)
		                {
		                  if(answer[e].ismatch && e!=q[0] && e!=q[1] && e!=q[2] && e!=q[3])
		                  {
		                    newin(0,0)=answer[e].in_match.col;
		                    newin(0,1)=answer[e].in_match.row;
		                    newin(0,2)=1;
		                    CImg<double> tempin=h*newin;
		                    tempin(0,0)/=tempin(0,2);
		                    tempin(0,1)/=tempin(0,2);
		                    tempin(0,2)/=tempin(0,2);
		                    dist=sqrt((tempin(0,0)-output_descriptors[i][e].col)*(tempin(0,0)-output_descriptors[i][e].col) + (tempin(0,1)-output_descriptors[i][e].row)*(tempin(0,1)-output_descriptors[i][e].row));
		                    no+=dist;
		                  }
		                }
		                if(best_no > no)
		                {
		                    best_h = h;
		                    best_no = no;
		                }
		            }
		            //best_h = best_h.get_invert();
		            
		            for(int e=0; e<output_descriptors[i].size(); e++)
		            { 
		              if(answer[e].ismatch)
		              {
		                  newin(0,0)=answer[e].in_match.col;
		                  newin(0,1)=answer[e].in_match.row;
		                  newin(0,2)=1;
		                  CImg<double> tempin=best_h*newin;
		                  tempin(0,0)/=tempin(0,2);
		                  tempin(0,1)/=tempin(0,2);
		                  tempin(0,2)/=tempin(0,2);
		                  dist=sqrt((tempin(0,0)-output_descriptors[i][e].col)*(tempin(0,0)-output_descriptors[i][e].col) + (tempin(0,1)-output_descriptors[i][e].row)*(tempin(0,1)-output_descriptors[i][e].row));
		                  if(dist>25)
		                  {  
		                      answer[e].ismatch = false;
		                      count_match--;
		                  }
		              }
		            }
		            feature_mark_append(input_image, output_images[i], appended_images[i], answer, output_descriptors[i]);
		            }
                
                //cout << "count_match: " << count_match << endl;
                
                if (argc > 3)
                {
                ratio = (double)count_match/(double)output_descriptors[i].size();
                for( j=desclist.size()-1; j>=0; j--)
		                if(count_match <= desclist[j].ct_match)
				                break;
		            j++;
		            temp.fname = argv[i+3];
		            temp.fname = temp.fname.substr( temp.fname.find_last_of("/")+1, string::npos);
		            temp.ct_match = count_match;
		            temp.ct_descriptors = output_descriptors[i].size();
		            temp.avg = ratio;
		            desclist.insert( desclist.begin()+j, temp);
		            out_name = out_path + inputFile + "_" + temp.fname;
		            //cout << out_name << endl;
		            appended_images[i].get_normalize(0, 255).save(out_name.c_str());
                }
	          }                
	          for(int i=0; i<desclist.size() && i<10; i++)
		            cout << desclist[i].fname << " " << desclist[i].ct_match << endl;
            cout << endl;
        }
        else if(part == "part3")
        {
            CImg<double> input_image(inputFile.c_str());
            string temporary = inputFile.substr( 0, inputFile.find_last_of('/'));
            char num = temporary[temporary.length()-2];
            inputFile = inputFile.substr( inputFile.find_last_of("/")+1, string::npos);
            inputFile = inputFile.substr( 0, inputFile.find('.'));
            cout << "input file: " << inputFile << endl;
            string out_name, out_path;
            
            if(argc==3)
            {
                out_path = "./output-a2-images/";
                out_name = out_path + inputFile + "-warped.png";
                warping(input_image,out_name);
                cout << out_name << endl;
            }
            else
            {
                out_path = "./output-a2-images/part3_images/seq1/";
                if(num == '2')
                    out_path = "./output-a2-images/part3_images/seq2/";
              
                CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
                vector < SiftDescriptor > descriptors = Sift::compute_sift(gray);
                
                kdtree *root;
                root = new kdtree;
                int level = 127;
                gen_kdtree( descriptors, &root, 0, level);
                
                vector < CImg<double> > output_images, output_gray;
                vector < vector < SiftDescriptor > > output_descriptors;
                for ( int of=3; of<argc; of++)
                {
                    CImg<double> temp(argv[of]);
                    output_images.push_back(temp);
                    output_gray.push_back(output_images[of-3].get_RGBtoHSI().get_channel(2));
                    output_descriptors.push_back(Sift::compute_sift(output_gray[of-3]));
                }
                
                int out_ct = output_descriptors.size(), count_match, sz, q[4];
                double thresh = 0.75, ratio;
                vector < ans > answer;
                
                for(int i=0; i<out_ct; i++)
                {
                    answer.clear();
                    count_match = check_out_image( &root, output_descriptors[i], thresh, answer);
                    //cout << "count_match: " << count_match << endl;
            
                    sz = output_descriptors[i].size();
                    
                    CImg<double> A(8,8), B(1,8), h(3,3), newin(1,3), best_h(3,3);
                    h(2,2)=1;
                    newin(0,2)=1;
                    double dist, no, best_no=FLT_MAX;
                    for(int u=0; u<500  && count_match >=4 ; u++)
                    {
                        do
                        {
                            q[0] = rand() % sz;
                        } while (!answer[q[0]].ismatch);
                        do
                        {
                            q[1] = rand() % sz;
                        }while(q[1] == q[0] || !answer[q[1]].ismatch);
                        do
                        {
                            q[2] = rand() % sz;
                        }while(q[2] == q[0] || q[2] == q[1] || !answer[q[2]].ismatch);
                        do
                        {
                            q[3] = rand() % sz;
                        }while(q[3] == q[0] || q[3] == q[1] || q[3] == q[2] || !answer[q[3]].ismatch);
                        for( int r=0; r<8; r++)
                        {
                            if(r%2 == 0)
                            {
                                A(0,r)=answer[q[r/2]].in_match.col;
                                A(1,r)=answer[q[r/2]].in_match.row;
                                A(2,r)=1;
                                for(int t=3; t<6; ++t)
                                    A(t,r)=0;
                                double x = -1 * output_descriptors[i][q[r/2]].col;
                                A(6,r)=(x * answer[q[r/2]].in_match.col);
                                A(7,r)=(x * answer[q[r/2]].in_match.row);
                                B(0, r) = output_descriptors[i][q[r/2]].col;
                                B(0, r+1) = output_descriptors[i][q[r/2]].row;
                            }
                            else
                            {   
                                for(int t=0; t<3; ++t)
                                    A(t,r)=0;
                                A(3,r)=answer[q[r/2]].in_match.col;
                                A(4,r)=answer[q[r/2]].in_match.row;
                                A(5,r)=1;
                                double y = -1 * output_descriptors[i][q[r/2]].row;
                                A(6,r)=(y * answer[q[r/2]].in_match.col);
                                A(7,r)=(y * answer[q[r/2]].in_match.row);
                            }
                        }
                        CImg<double> tA = A.invert();
                        CImg<double> C=tA*B;
                        for(int tr=0; tr<8; tr++)
                            h(tr%3,tr/3)=C(0,tr);
                        no=0;
                        newin(0,2)=1;
                        for(int e=0; e<sz; e++)
                        {
                            if(answer[e].ismatch && e!=q[0] && e!=q[1] && e!=q[2] && e!=q[3])
                            {
                                newin(0,0)=answer[e].in_match.col;
                                newin(0,1)=answer[e].in_match.row;
                                CImg<double> tempin=h*newin;
                                tempin(0,0)/=tempin(0,2);
                                tempin(0,1)/=tempin(0,2);
                                tempin(0,2)/=tempin(0,2);
                                dist=sqrt((tempin(0,0)-output_descriptors[i][e].col)*(tempin(0,0)-output_descriptors[i][e].col) + (tempin(0,1)-output_descriptors[i][e].row)*(tempin(0,1)-output_descriptors[i][e].row));
                                no+=dist;
                            }
                        }
                        if(best_no > no)
                        {
                            best_h = h;
                            best_no = no;
                        }
                    }
                    best_h = best_h.get_invert();
                    int a,b, r= output_images[i]._height, c= output_images[i]._width;
                    CImg<double> tempout = output_images[i];
                    newin(0,2)=1;
                    for(int tr=0; tr<r; tr++)
                        for(int tc=0; tc<c; tc++)
                        {
                            tempout(tc,tr,0,0)=0;
                            tempout(tc,tr,0,1)=0;
                            tempout(tc,tr,0,2)=0;
                            newin(0,0)=tc;
                            newin(0,1)=tr;
                            newin=best_h*newin;
                            newin(0,0)/=newin(0,2);
                            newin(0,1)/=newin(0,2);
                            newin(0,2)/=newin(0,2);
                            if(newin(0,0) >= 0 && newin(0,1) >= 0 && newin(0,1)<output_images[i]._height && newin(0,0)<output_images[i]._width)
                            {
                                tempout(tc,tr,0,0)=output_images[i](newin(0,0),newin(0,1),0,0);
                                tempout(tc,tr,0,1)=output_images[i](newin(0,0),newin(0,1),0,1);
                                tempout(tc,tr,0,2)=output_images[i](newin(0,0),newin(0,1),0,2);                    
                            }
                        }
                     
                    string tempfname = argv[i+3];
                    tempfname = tempfname.substr( tempfname.find_last_of("/")+1, string::npos);
                    tempfname = tempfname.substr( 0, tempfname.find('.'));
                    
                    out_name = out_path + tempfname + "-warped.png";
                    tempout.get_normalize(0, 255).save(out_name.c_str());
                    cout << out_name << endl;
                    
                }
            }
        }
        else
            throw std::string("unknown part!");
    
    }
    catch(const string &err) 
    {
        cerr << "Error: " << err << endl;
    }
}

/*
 for(int i=0; i<descriptors.size(); i++)
 {
 //cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
 //for(int l=0; l<128; l++)
 //    cout << descriptors[i].descriptor[l] << "," ;
 //cout << ")" << endl;
 
 for(int j=-1; j<=1; j++)
 for(int k=-1; k<=1; k++)
 if(j==0 || k==0)
 for(int p=0; p<3; p++)
 {    
 if(descriptors[i].col+k < input_image.width() && descriptors[i].row+j < input_image.height() && descriptors[i].col+k >= 0 && descriptors[i].row+j >= 0)
 input_image(descriptors[i].col+k, descriptors[i].row+j, 0, p)=0;
 }
 }
*/

void gen_kdtree(vector < SiftDescriptor > descriptors, kdtree **root, int i, int level)
{
    //cout << i << endl;
    if(i < level)
    {
        int s = descriptors.size();
        int mid = s/2, left, right, index;
        std::sort (descriptors.begin(), descriptors.end(), compare(i));
        index = -1;
        left = mid;
        right = mid;
        while(left>=0 && right<s)
        {
            if(descriptors[left].descriptor[i] != descriptors[mid].descriptor[i])
            {   
                index = left+1;
                break;
            }
            else if(descriptors[right].descriptor[i] != descriptors[mid].descriptor[i])
            {
                index = right;
                break;
            }
            if(descriptors[left].descriptor[i] == descriptors[mid].descriptor[i])
                left--;
            if(descriptors[right].descriptor[i] == descriptors[mid].descriptor[i])
                right++;
        }
        if(index == -1)
            index = mid;
    
        (*root)->val = descriptors[index].descriptor[i];
        (*root)->index = i;
        vector < SiftDescriptor > ::iterator middle = descriptors.begin() + index;
        vector < SiftDescriptor > ld( descriptors.begin(), middle);
        vector < SiftDescriptor > rd( middle, descriptors.end());
        if(ld.size()+rd.size() == s)
        {
            //cout << " left_size: " << ld.size() << " right_size: " << rd.size() << " val: " << root->val << " index: " << root->index << " dsize: " << root->descriptors.size() << endl;
            if(ld.size() == 0 || rd.size() == 0)
            {
                (*root)->val = -1;
                (*root)->index = -1;
                (*root)->l = NULL;
                (*root)->r = NULL;
                (*root)->descriptors = descriptors;
                //cout << i<< " currsize: " << descriptors.size() << " dsize: " << root->descriptors.size() << endl;
            }
            else
            {
                (*root)->l = new kdtree;
                gen_kdtree( ld, &((*root)->l), i+1, level);
                (*root)->r = new kdtree;
                gen_kdtree( rd, &((*root)->r), i+1, level);
            }
        }
    }
    else
    {
        (*root)->val = -1;
        (*root)->index = -1;
        (*root)->l = NULL;
        (*root)->r = NULL;
        (*root)->descriptors = descriptors;
        //cout << i<< " currsize: " << descriptors.size() << " dsize: " << root->descriptors.size() << endl;
    }
} 

void print_kdtree(kdtree **root)
{
    if(*root != NULL)
    {
        if((*root)->index == -1)
        {
            int s = (*root)->descriptors.size();
            for(int i=0; i<s; i++)
                cout << "Descriptor #" << i << ": " << (*root)->descriptors[i].descriptor.size() << endl;
        }
        else
        {
            cout << (*root)->val << " " << (*root)->index << " " << (*root)->descriptors.size() << endl;
            print_kdtree( &((*root)->l));
            print_kdtree( &((*root)->r));
        }
    }
}

double euclidean_dist(SiftDescriptor &a, SiftDescriptor &b)
{
    double ans = 0;
    for(int i=0; i<128; i++)
        ans+= ((a.descriptor[i]-b.descriptor[i]) * (a.descriptor[i]-b.descriptor[i]));
    ans = sqrt(ans);
    return ans;
}

void search_kdtree(kdtree **root, SiftDescriptor curr, double &min, double temp, SiftDescriptor &match)
{
    if(*root != NULL)
    {
        if((*root)->index == -1)
        {
            int s = (*root)->descriptors.size();
            for(int i=0; i<s; i++)
            {
                double ed = euclidean_dist (curr, (*root)->descriptors[i]);
                if (ed < min && ed > temp)
                {
                    min = ed;
                    match = (*root)->descriptors[i];
                }
            }
        }
        else
        {
            int i = (*root)->index;
            //if(curr.descriptor[i] < (*root)->val)
            search_kdtree( &((*root)->l), curr, min, temp, match);
            //else
            search_kdtree( &((*root)->r), curr, min, temp, match);
        }
    }
}

int check_out_image(kdtree **root, vector < SiftDescriptor > out_descriptors, double thresh, vector < ans > &answer)
{
    int s = out_descriptors.size(), count_match = 0, i, j;
    double min, min2;
    SiftDescriptor close, close2;
    ans temp;
    for(i=0; i<s; i++)
    {
        min2 = -1;
        min = FLT_MAX;
        search_kdtree( root, out_descriptors[i], min, min2, close);
        min2 = FLT_MAX;
        search_kdtree( root, out_descriptors[i], min2, min, close2);
        //cout << min << " " << min2 << " " << min/min2 << endl;
        temp.in_match = close;
        temp.mined = min;
        temp.min2ed = min2;
        temp.val = (min/min2);
        temp.ismatch = false;
        if( temp.val <= thresh )
            temp.ismatch = true;
        if(temp.ismatch)
            count_match++;
        //cout << temp.mined << " " << temp.min2ed << " " << temp.val << " " << temp.ismatch << endl;
        for(j=0; j<answer.size() && temp.ismatch; j++)
            if(answer[j].in_match.row == temp.in_match.row && answer[j].in_match.col == temp.in_match.col && answer[j].ismatch)
                break;
        if(j < answer.size() && temp.ismatch && answer[j].ismatch)
        {
            count_match--;
            if(answer[j].mined > temp.mined)
                answer[j].ismatch = false;
            else
                temp.ismatch = false;
        }
        answer.push_back(temp);
    }
    return count_match;
}
  
/*
void feature_mark(CImg<double> &input_image, vector < SiftDescriptor > &descriptors)
{
    for(int i=0; i<descriptors.size(); i++)
    {
        for(int j=-1; j<=1; j++)
            for(int k=-1; k<=1; k++)
                if(j==0 || k==0)
                    for(int p=0; p<3; p++)
                    {    
                        if(descriptors[i].col+k < input_image.width() && descriptors[i].row+j < input_image.height() && descriptors[i].col+k >= 0 && descriptors[i].row+j >= 0)
                            input_image(descriptors[i].col+k, descriptors[i].row+j, 0, p)=0;
                    }
    }
}
*/

void feature_mark_append(CImg<double> input_image, CImg<double> &output_image, CImg<double> &appended_image, vector <ans> &answer, vector < SiftDescriptor > &output_descriptors)
{
    int s = output_descriptors.size();
    for(int i=0; i<s; i++)
        if(answer[i].ismatch)
            for(int j=-1; j<=1; j++)
                for(int k=-1; k<=1; k++)
                    if(j==0 || k==0)
                        for(int p=0; p<3; p++)
                        {    
                            if(output_descriptors[i].col+k < output_image.width() && output_descriptors[i].row+j < output_image.height() && output_descriptors[i].col+k >= 0 && output_descriptors[i].row+j >= 0)
                                output_image(output_descriptors[i].col+k, output_descriptors[i].row+j, 0, p)=0;
                            if(answer[i].in_match.col+k < input_image.width() && answer[i].in_match.row+j < input_image.height() && answer[i].in_match.col+k >= 0 && answer[i].in_match.row+j >= 0)
                                input_image(answer[i].in_match.col+k, answer[i].in_match.row+j, 0, p)=0;
                        }
    int firstrow = input_image._height, firstcol = input_image._width;
    appended_image = input_image;
    appended_image.append( output_image, 'x');
    const unsigned char color[] = { 255, 255, 0};
    for(int i=0; i<s; i++)
        if(answer[i].ismatch)
            appended_image.draw_line(answer[i].in_match.col, answer[i].in_match.row, output_descriptors[i].col + firstcol, output_descriptors[i].row, color);
}

CImg<double> matrixmult( CImg<double> A, CImg<double> B)
{
    CImg<double> C(1,8);
    //cout << endl;
    double val;
    for(int ii=0;ii<8;ii++)
    {   
        C(0,ii)=0.0;
        for(int jj=0;jj<8;jj++)
        {
            //cout << C(0,ii) << endl;
            //cout << A(jj,ii) << " " << B(0,jj) << " " << (A(jj,ii)*B(0,jj)) << endl;
            val = (A(jj,ii)*B(0,jj));
            C(0,ii) = C(0,ii) + val;
        }
        //cout << C(0,ii) << endl;
    }
    //cout << endl;
    //for(int ij = 0; ij < 8; ij++)
    //    cout << C(0,ij) << " ";
    //cout << endl;
    //int d;
    //cin>>d;
    return C;
}
void warping(CImg<double> input_image, string out_name)
{
    CImg<double> h(3,3), inp(1,3);
    h(0,0)=0.907;
    h(1,0)=0.258;
    h(2,0)=-182;
    h(0,1)=-0.153;
    h(1,1)=1.44;
    h(2,1)=58;
    h(0,2)=-0.000306;
    h(1,2)=0.000731;
    h(2,2)=1;
    
    h = h.get_invert();
    
    inp(0,2)=1;
    int r = input_image._height, c = input_image._width, a, b;
    CImg<double> output_image = input_image;
    for(int tr=0; tr<r; tr++)
        for(int tc=0; tc<c; tc++)
        {
            output_image(tc,tr,0,0)=0;
            output_image(tc,tr,0,1)=0;
            output_image(tc,tr,0,2)=0;
            inp(0,0)=tc;
            inp(0,1)=tr;
            inp=h*inp;
            //cout << inp(0,0) << " " << inp(0,1) << " " << inp(0,2) << endl;
            inp(0,0)/=inp(0,2);
            inp(0,1)/=inp(0,2);
            inp(0,2)/=inp(0,2);
            //cout << inp(0,0) << " " << inp(0,1) << " " << inp(0,2) << endl;
            if(inp(0,0) >= 0 && inp(0,1) >= 0 && inp(0,1)<r && inp(0,0)<c)
            {
                    //cout << a << " " << b << " " << input_image(a,b,0,p);
                    output_image(tc,tr,0,0)=input_image(inp(0,0),inp(0,1),0,0);
                    output_image(tc,tr,0,1)=input_image(inp(0,0),inp(0,1),0,1);
                    output_image(tc,tr,0,2)=input_image(inp(0,0),inp(0,1),0,2);                    
                    //cout << " " << output_image(tc,tr,0,p) << endl;
            }
        }
        
    //cout << output_image(909,605,0,2) << endl;
    output_image.get_normalize(0, 255).save(out_name.c_str());
}