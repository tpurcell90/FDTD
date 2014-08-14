#include "toBitMap.hpp"

#include <stdio.h>
#include <algorithm>
#include <complex>

using namespace std;

void fieldToBitMap (FDTDField &F, string filename)
{
  FILE *outMap;
  int w = F.nx();
  int h = F.ny();
  int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

  vector<char> img(3*w*h);

  int r,g,b;

  for(int i = 0; i < w; i++)
  {
    for(int j = 0; j < h; j++)
    {
      int x = i;
      int y = (h-1) - j;
      double max = real(*max_element(F.Ez_->data(),F.Ez_->data()+h*w,[](complex<double> ii, complex<double> jj){return real(ii) < real(jj);}));
      double min = real(*min_element(F.Ez_->data(),F.Ez_->data()+h*w,[](complex<double> ii, complex<double> jj){return real(ii) < real(jj);}));

      //Adjust this to tune the color
      r = int(255.0*real((min + F.Ez_->point(i,j)))/(max-min));// //red[i][j]*255;
      g = int(255.0*real((min + F.Ez_->point(i,j)))/(max-min)); ;//green[i][j]*255;
      b = 0;//blue[i][j]*255;
      if (r > 255) r=255;
      // if (g > 255) g=255;
      // if (b > 255) b=255;

      img[(x+y*w)*3+2] = (unsigned char)(r);
      img[(x+y*w)*3+1] = (unsigned char)(g);
      img[(x+y*w)*3+0] = (unsigned char)(b);
    }
  }
  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);

  bmpinfoheader[ 4] = (unsigned char)(w);
  bmpinfoheader[ 5] = (unsigned char)(w >> 8);
  bmpinfoheader[ 6] = (unsigned char)(w >> 16);
  bmpinfoheader[ 7] = (unsigned char)(w >> 24);
  bmpinfoheader[ 8] = (unsigned char)(h);
  bmpinfoheader[ 9] = (unsigned char)(h >> 8);
  bmpinfoheader[10] = (unsigned char)(h >> 16);
  bmpinfoheader[11] = (unsigned char)(h >> 24);

  outMap = fopen((char*)filename.c_str(),"wb");
  fwrite(bmpfileheader,1,14,outMap);
  fwrite(bmpinfoheader,1,40,outMap);
  for(int i=0; i<h; i++)
  {
      fwrite(img.data()+(w*(h-i-1)*3),3,w,outMap);
      fwrite(bmppad,1,(4-(w*3)%4)%4,outMap);
  }
  fclose(outMap);
}
