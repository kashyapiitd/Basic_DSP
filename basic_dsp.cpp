#include <iostream>
#include <vector>

using namespace std;

/*
Performs 1-D Convolution between two Data Vectors 
*/

void perform_convolution(vector<int> & x, vector<int> & h)
{
    vector<int> y (x.size() + h.size() - 1, 0);
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
    // Print the Convolved Output
    for (int i = 0 ; i < y.size(); i++)
    {
        cout << y[i] << " ";
    }
    cout << endl;
}
int main()
{
	vector<int> x = {1,2,3,4,5};
	vector<int> h = {1,2,1,2,1};
	perform_convolution(x,h);
	return 1;
}
