#include <iostream>
#include <unistd.h>

using namespace std;

int main()
{
    
    cout << endl;
    cout << "  test progress bar v1.0" << endl;
    
    for(int i = 0; i < 100; ++i)
    {   
        if(i < 10) { cout << "  [00" << i << "] "; }
        
        else if(i < 100) { cout << "  [0" << i << "] "; }
        
        else { cout << "  [" << i << "] "; }
        
        for(int j = 0; j < i; ++j)
        {    cout << ">";    }
        
        cout.flush();
        
        usleep(50000);
        
        cout << "\r";
    }
    
    cout << "\r";
    cout << "  [100] Done                     ";
    
    cout << endl;
    cout << endl;
    
    return 0;
}


