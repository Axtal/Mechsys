/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
 * Copyright (C) 2007 Daichao Sheng, David A. F. Oliveira               *
 * Copyright (C) 2007 Kristian Krabbenhøft                              *
 * Copyright (C) 2007 Marcos A. A. Santos, Jidong Zhao                  *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/* Matrix-Vector expression template library check
 *
 */

// STL
#include <iostream>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>

using namespace std;
#ifndef USE_MTL4
  using namespace LinAlg;
#endif

int main(int argc, char **argv) try
{
#ifdef USE_MTL4
    cout << "This test only works with Raul's expression template library. Not with MTL4\n";
#else
    cout << "Matrix-Vector expression template library check" << endl;
    cout << endl;
    bool chk_success = true;

    cout << "Vector addition... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A = 1, 1, 1;
        Vector<double> B(3); B = 1, 1, 1;
        Vector<double> R(3);
        R = A + B + A;
        if (R(0) != 3) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Vector subtraction... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A.SetValues(1);
        Vector<double> B(3); B.SetValues(1);
        Vector<double> R(3);
        R = A - B;
        if (R(0) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else        cout << "not passsed. " << endl << endl;
    }

    cout << "Vector multiplication by scalar... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A.SetValues(1);
        Vector<double> B(3); B.SetValues(1);
        Vector<double> R(3);
        R = 2*A + B;
        if (R(0) != 3) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else        cout << "not passsed. " << endl << endl;
    }

    cout << "Vector division by scalar... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A.SetValues(2);
        Vector<double> B(3); B.SetValues(1);
        Vector<double> R(3);
        R = A/2 + B;
        if (R(0) != 2) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else        cout << "not passsed. " << endl << endl;
    }

    cout << "Vector tranpose... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A = 0, 1, 2;
        Matrix<double> R(3,1);
        R = trn(A);
        if (R(0,0) != 0 || R(0,1) != 1 || R(0,2) != 2) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix addition... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); A.SetValues(1);
        Matrix<double> B(2,3); B.SetValues(1);
        Matrix<double> R;
        R = A + B + A;
        if (R(0,0) != 3) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix subtraction... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); A.SetValues(1);
        Matrix<double> B(2,3); B.SetValues(1);
        Matrix<double> R;
        R = A - B;
        if (R(0,0) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix multiplication by scalar... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); A.SetValues(1);
        Matrix<double> B(2,3); B.SetValues(1);
        Matrix<double> R;
        R = 2*A + B;
        if (R(0,0) != 3 || R(1,2) != 3) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix division by scalar... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); A.SetValues(3);
        Matrix<double> B(2,3); B.SetValues(1);
        Matrix<double> R;
        R = A/3 + B;
        if (R(0,0) != 2 || R(1,2) != 2) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix tranpose... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); 
        A = 0, 1, 2,
            3, 4, 5;
        Matrix<double> R;
        R = trn(A);
        if (R(0,0) != 0 || R(1,0) != 1 || R(2,0) != 2) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Vector Matrix multiplication... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); 
        Matrix<double> B(1,3); 
        A = 1, 1, 1;
        B = 0, 1, 2;
        Matrix<double> R;
        R = A*B;
        if (R(0,0) != 0 || R(0,1) != 1 || R(0,2) != 2) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix Vector multiplication... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); 
        Vector<double> B(3); 
        A = 0, 1, 2,
            3, 4, 5;
        B = 1, 1, 1;
        Vector<double> R;
        R = A*B;
        //std::cout << "Preview = " << A << std::endl;
        //std::cout << "Result  = " << R << std::endl;
        if (R(0) != 3 || R(1) != 12) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix Matrix multiplication... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); 
        Matrix<double> B(3,2); 
        A = 0, 1, 2,
            3, 4, 5;
        B = 1, 1,
            1, 1,
            1, 1;
        Matrix<double> R;
        R = A*B;
        if (R(0,0) != 3 || R(1,0) != 12) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix determinant... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,3); 
        A = 2, 0, 0,
            0, 2, 0,
            0, 0, 2;
        double R;
        R = Det(A);
        if (R != 8) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }
    
    cout << "Matrix inversion... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(4,4); 
        A = 2, 0, 0, 0,
            0, 2, 0, 0,
            0, 0, 2, 0,
            0, 0, 0, 2;
        Matrix<double> R;
        Inv(A, R);
        if (R(0,0) != 0.5 || R(1,1) != 0.5 || R(2,2) != 0.5) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix by transposed Matrix multiplication... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(2,3); 
        Matrix<double> B(2,3); 
        A = 0, 1, 2,
            3, 4, 5;
        B = 1, 1, 1,
            1, 1, 1;
        Matrix<double> R;
        R = A*trn(B);
        if (R(0,0) != 3 || R(1,0) != 12) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "transposed Matrix by Matrix multiplication... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,2); 
        Matrix<double> B(3,2); 
        A = 0, 3,
            1, 4,
            2, 5;
        B = 1, 1,
            1, 1,
            1, 1;
        Matrix<double> R;
        R = trn(A)*B;
        if (R(0,0) != 3 || R(1,0) != 12) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "transposed Matrix by transposed Matrix multiplication... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,2); 
        Matrix<double> B(2,3); 
        A = 0, 3,
            1, 4,
            2, 5;
        B = 1, 1, 1,
            1, 1, 1;
        Matrix<double> R;
        R = trn(A)*trn(B);
        //std::cout << "Result  = " << R << std::endl;
        if (R(0,0) != 3 || R(1,0) != 12) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Vector by transposed Vector multiplication... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3); A.SetValues(2);
        Vector<double> B(3); B.SetValues(2);
        Matrix<double> R;
        R = A*trn(B);
        if (R(0,0) != 4 || R(2,2) != 4) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Vector expression ( R = A + 2*B - (A + B + B) ) evaluation... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3);
        Vector<double> B(3);
        A = 1, 2, 3;
        B = 4, 5, 6;
        Vector<double> R;
        R = A + 2*B - (A + B + B);
        if (R(0) != 0 || R(1) != 0 || R(2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix expression ( R = 2*A - A*B*3*2 - (A + A - 2*3*trn(trn(A) * trn(B)) ) + trn(A)*B - trn(trn(B) * A) ) evaluation... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,3); 
        Matrix<double> B(3,3); 
        A = 0, 1, 2,
            3, 4, 5,
            6, 7, 8;
        B = A;
        Matrix<double> R;
        R = 2*A - A*B*3*2 - (A + A - 2*(3*trn(trn(A) * trn(B))) ) + trn(A)*B - trn(trn(B) * A);
        if (R(0,0) != 0 || R(1,0) != 0 || R(2,2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Vector expression ( R -= 2*trn(A) - trn(A + A) ) evaluation... " << endl;
    {
        bool         passed = true;
        Vector<double> A(3);
        A(0) = 1; A(1) = 2; A(2) = 3; 
        Matrix<double> R(1,3);
        R.SetValues(0);
        R -= 2*trn(A) - trn(A + A);
        if (R(0,0) != 0 || R(0,1) != 0 || R(0,2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix expression ( R += 2*A*B*trn(A) - trn(2*A*trn(B*A)) ) evaluation... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,3); 
        Matrix<double> B(3,3); 
        A = 0, 1, 2, 
            3, 4, 5, 
            6, 7, 8; 
        B = A;
        Matrix<double> R(3,3);
        R.SetValues(0);
        R += 2*A*B*trn(A) - trn(2*A*trn(A*B));
        if (R(0,0) != 0 || R(1,0) != 0 || R(2,2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix expression ( R = Det(A)*trn(inv(B)) - trn(Det(A)*inv(B)) ) evaluation... " << endl;
    {
        bool         passed = true;
        Matrix<double> A(3,3); 
        Matrix<double> B(3,3), Bi; 
        A = 2, 2, 3, 
            4, 5, 6, 
            7, 8, 9; 
        B = A;
        Inv (B, Bi);
        Matrix<double> R(3,3);
        R = Det(A)*trn(Bi) - trn(Det(A)*Bi);
        //std::cout << "Result  = " << R << std::endl;
        if (R(0,0) != 0 || R(1,0) != 0 || R(2,2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix expression ( R = trn(A)*B*trn(C) - trn(C*trn(B)*A) evaluation... " << endl;
    {                                      
        bool         passed = true;
        Matrix<double> A(2,3); 
        Vector<double> B(2);
        Vector<double> C(3);
        A = 1, 2, 3, 
            4, 5, 6; 
        B.SetValues(2);
        C.SetValues(3);
        Matrix<double> R;
        R = trn(A)*B*trn(C) - trn(C*trn(B)*A);
        if (R(0,0) != 0 || R(1,0) != 0 || R(2,2) != 0) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    cout << "Matrix expression ( R += trn(B)*D*B*Det(J)/c  - trn(trn(D*B)*B)*Det(J)/c evaluation... " << endl;
    {                                      
        bool         passed = true;
        Matrix<double> B(2,3); 
        Matrix<double> D(2,2); 
        Matrix<double> J(3,3); 
        B = 1, 2, 3, 
            4, 5, 6; 
        D = 2, 3, 
            4, 5; 
        J = 1, 0, 0,
            0, 1, 0,
            0, 0, 1;
        double c = 1;
        Matrix<double> R(3,3);
        R.SetValues(1);
        R+= trn(B)*D*B*Det(J)/c - trn(trn(D*B)*B)*Det(J)/c;
        if (R(0,0) != 1 || R(1,0) != 1 || R(2,2) != 1) passed = false;
        if (passed) cout << "ok. " << endl << endl;
        else
        {   
            chk_success = false;
            cout << "not passsed. " << endl << endl; 
        }
    }

    if (chk_success) cout << "Check result: ok. "          << endl << endl;
    else             cout << "Check result: not passsed. " << endl << endl;

#endif
    return 0;
}
MECHSYS_CATCH
