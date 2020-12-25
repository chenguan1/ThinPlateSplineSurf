using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;


namespace TPSSurfTest
{
    class Program
    {
        [DllImport(@"E:\CodeProjects\C++\TPSDemo\x64\Release\TPSSurf.dll", CharSet = CharSet.Auto, CallingConvention = CallingConvention.Cdecl)]
        private static extern int CalTpsSurf(int ptcount, double[] xs, double[] ys, double[] zs, int rows, int cols, double[] gridOut);

        
        static void Main(string[] args)
        {
            double[] xs = new double[] { 0, 99, 0 };
            double[] ys = new double[] { 0, 0, 49 };
            double[] zs = new double[] { 0, 100, 50 };

            double[] grid = new double[100 * 50];
            int rtn = CalTpsSurf(3, xs, ys, zs, 50, 100, grid);

            Console.WriteLine(rtn);

        }
    }
}
