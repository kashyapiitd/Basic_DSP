/*
License: MIT License (http://www.opensource.org/licenses/mit-license.php)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
/* (C) 2018 Dhritiman Kashyap */

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "log.h"

using namespace std;

/*

Function: perform_linear_convolution
Input: input signal x and filter coefficients h.
Output: output convolved signal y.
Description: Performs 1-D Linear Convolution between two Data Vectors 

*/

template<typename T>
vector<T> perform_linear_convolution(vector<T> & x, vector<T> & h)
{
    vector<T> y (x.size() + h.size() - 1, 0);
    for (int i = 0; i < y.size(); i++)
    {
        for (int j = 0; i-j >= 0 ; j++)
        {
            if ( j >= x.size() || (i-j >= h.size()))
            {
                continue;
            }
            y[i] += x[j] * h[i-j];
        }
    }
    return y;
}


/*

Function: perform_circular_convolution
Input: input signal x and filter coefficients h.
Output: output convolved signal y.
Description: Performs 1-D Circular Convolution between two Data Vectors 

*/
template<typename T>
vector<T> perform_circular_convolution (vector<T> x, vector<T> h)
{
    int output_length = max(x.size(), h.size());
    if (x.size() != output_length)
    {
        x.resize(output_length,0);
    }
    if (h.size() != output_length)
    {
        h.resize(output_length, 0);
    }
    vector<T> y (output_length, 0);
    for (int n = 0; n < output_length; n++)
    {
        for (int k = 0; k < output_length; k++)
	{
	    if (n < k )
	    {
	        y[n] += x[k] * h[n+output_length-k];
	    }
	    else
	    {
	        y[n] += x[k] * h[n-k];
	    }
	}
    }
    return y;
}

/*
Function: overlap_and_save_filtering.
Input: input signal x and filter coefficients h.
Output: output signal y.
Description: Following are the stepsinvolved in overlap and save filtering
    1. Divide the incoming signal into blocks of lenght L (where L > M, M being the size of the filter coefficients) and prepend the block with M-1 samples of the previous block, thus making the 
       size of each block as L+M-1. For the first block , prepend M-1 zeros.
    2. Perform Circular convolution between a block x_l[n] and the filter coefficients h[n]. The result of the circular convolution is y_l[n] with length max(L+M-1,M) which is L+M-1.
    3. Ignore the first M-1 samples and accumulate the rest L samples for each block to get the final result.
This approach is used generally when we can employ circular convolution.
*/

template<typename T>
vector<T> overlap_and_save_filtering(vector<T>& x, vector<T>& h, int block_size)
{
    int M = h.size();
    int L = block_size;
    int num_blocks = x.size() / L;
    if (x.size() % L != 0)
    {
        num_blocks++;
    }
    int output_length = x.size() + h.size() - 1;
    vector<T> y(output_length, 0);
    vector<T> input_block (L+M-1, 0);
    vector<T> output_block (L+M-1, 0);
    vector<T> overlap_buffer (M-1, 0);
    int start_index = 0;
    int end_index = L;
    while(num_blocks > 0)
    {
        copy(overlap_buffer.begin(), overlap_buffer.end(), input_block.begin());
	copy(x.begin() + start_index, x.begin() + end_index, input_block.begin() + M - 1);
	// Perform Circular convolution
        output_block = perform_circular_convolution<T> (input_block, h);
	copy(output_block.begin() + M - 1, output_block.end(), y.begin() + start_index);
	// Fill the overlap buffer for the next iteration
	copy(input_block.begin() + L, input_block.end(), overlap_buffer.begin());
	start_index += L;
	end_index += L;
	num_blocks--;
    }
    return y;



}

/*
Function: overlap_and_add_filtering.
Input: input signal x and filter coefficients h (of length M).
Output: output signal y.
Description: Following are the steps involved in overlap and add filtering
    1. Divide the signal into blocks of length L with no overlap. The first block starts at index 0 and extends till L-1, second block starts at L till 2L-1 and so on.
    2. Perform Linear convolution between each input block x_l[n] and the filter coefficients h[n] resulting into an output block of length L+M-1.
    3. Accumulate these output blocks in such a way that each output block overlaps with the previous output block by M-1 samples. To get the resulting sample for the overlap region, add the respective 
       overlapping samples.
This approach is used generally when we can employ linear convolution.
*/
template<typename T>
vector<T> overlap_and_add_filtering(vector<T>& x, vector<T>& h, int block_size)
{
    int M = h.size(); // M is the Filter Length
    int L = block_size; // L is the Input Signal Block Length
    int num_blocks = x.size() / L;
    if (x.size() % L != 0)
    {
        num_blocks++;
    }
    int output_length = x.size() + h.size() -1;
    LOG("Final Output Buffer Size: ", output_length);
    vector<T> y(output_length, 0); // Initialise the output as a vector of zeros
    vector<T> input_block;
    vector<T> output_block;
    int start_index = 0;
    int end_index =  L;
    while(num_blocks > 0)
    {
	// Resize the input_block vector to be a vector of zeros and of size L
	input_block.resize(L,0);
        // Copy one input block
        copy(x.begin() + start_index, x.begin() + end_index, input_block.begin());
	// Perform Filtering of the Block
	output_block = perform_linear_convolution<T>(input_block, h);
	// Write result back in y
	transform(output_block.begin(), output_block.end(), y.begin() + start_index, y.begin() + start_index, [](int x, int y){return x + y;});
	// Clear input_block fornext iteration
	input_block.clear();
	// Update start_index and end_index
	start_index += L;
	end_index += L;
	num_blocks--;
    }
    return y;	
}

/*
Function: perform_real_time_filtering.
input: input sample x and the filter coeff.
output: filtererd output sample.
Description: The function simulates a real time filtering scenario in which filtering is done on a sample by sample basis as opposed to a block by block basis as done in 
             overlap add or overlap save method. We internally maintain a data buffer which is updated with every new input sample and filtering operation is performed.

*/

template<typename T>
T perform_real_time_filtering(T sample, vector<T> &filter_coeff)
{
    // Initialise the data_buffer to be a zero buffer with size same as the size of filter coefficient
    static vector<T> data_buffer(filter_coeff.size(), 0);
    // Always insert the incoming sample into the 0th position.
    data_buffer[0] = sample;
    for (int i = 0 ; i < data_buffer.size(); i++)
    {
        cout << data_buffer[i] << " ";
    }
    cout << endl;
    // Perform the Actual filtering
    T output_sample = 0;
    for (int i = 0; i < filter_coeff.size(); i++)
    {
        output_sample += filter_coeff[i] * data_buffer[i];
    }
    // Update the data_buffer
    for (int i= data_buffer.size() - 2; i >= 0; i--)
    {
        data_buffer[i+1] = data_buffer[i];
    }
    return output_sample;

} 

/*
Function: Main
Description: Driver Function simulating all different functions.

*/
int main()
{
	vector<int> x = {1,2,3,2};
	vector<int> h = {1,2,1};
	//vector<int> y = perform_linear_convolution<int>(x,h);
	//vector<int> y = overlap_and_add_filtering<int>(x,h,3);
	//vector<int> y = perform_circular_convolution<int>(x,h);
	//vector<int> y = overlap_and_save_filtering<int>(x,h,3);
        // Perform Real Time Filtering
	vector<int> y;
        for (int i = 0 ; i < x.size(); i++)
	{
            y.push_back(perform_real_time_filtering<int>(x[i], h));
	}
        // Print the Convolved Output
        for (int i = 0 ; i < y.size(); i++)
        {
            cout << y[i] << " ";
        }
        cout << endl;
	return 1;
}
