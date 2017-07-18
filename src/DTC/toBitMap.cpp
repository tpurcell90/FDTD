#include "toBitMap.hpp"

int toGValue(double a)
{
    int a_ind = int( std::floor( std::abs(a) * 255 + 0.5) ) ;
    if(a_ind >= 255)
        return int( 255 * G_vals[255]);
    else
        return int( 255 * G_vals[a_ind]);
}

int toRValue(double a)
{
    int a_ind = int( std::floor( std::abs(a) * 255 + 0.5) ) ;
    if(a_ind >= 255)
        return int( 255 * R_vals[255]);
    else
        return int( 255 * R_vals[a_ind]);
}
int toBValue(double a)
{
    int a_ind = int( std::floor( std::abs(a) * 255 + 0.5) ) ;
    if(a_ind >= 255)
        return int( 255 * B_vals[255]);
    else
        return int( 255 * B_vals[a_ind]);
}



void GridToBitMap (std::vector<int_grid_ptr> o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);

    int min = 0, max = 0;
    int temp1 = 0, temp2 = 0;
    if(part == PLOTTYPE::POW)
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxReal(grid);
            min += temp1*temp1;
            max += temp2*temp2;
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxReal(grid);
            min += temp1*temp1;
            max += temp2*temp2;
        }
        max = log(max);
        min = log(min);
    }
    else
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxReal(grid);
            min += temp1;
            max += temp2;
        }
    }
    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    double diff = max - min;
    if(part == PLOTTYPE::POW)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += pow(grid->point(i+loc[0],j+loc[1]), 2.0);
                val = (val - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += pow(grid->point(i+loc[0],j+loc[1]), 2.0);
                val = (log(val) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
                int x = 0, y = 0;
        double val = 0.0;
        for(int ii = 0; ii < w; ii++)
        {
            for(int jj = 0; jj < h; jj++)
            {
                val = 0;
                x = ii; y = (h-1) - jj;
                for(auto & grid : o)
                    val += (grid->point(ii+loc[0],jj+loc[1]) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);

    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (std::vector<real_grid_ptr> o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);

    double min   = 0,     max = 0;
    double temp1 = 0.0, temp2 = 0.0;
    if(part == PLOTTYPE::POW)
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxPower(grid);
            min += temp1*temp1;
            max += temp2*temp2;
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxPower(grid);
            min += temp1*temp1;
            max += temp2*temp2;
        }
        min = log(min);
        max = log(max);
    }
    else
    {
        for(auto & grid: o)
        {
            std::tie(temp1, temp2) = findMinMaxReal(grid);
            min += temp1;
            max += temp2;
        }
    }
    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    double diff = max - min;
    if(part == PLOTTYPE::POW)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += pow(grid->point(i+loc[0],j+loc[1]), 2.0);
                val = (val - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += pow(grid->point(i+loc[0],j+loc[1]), 2.0);
                val = ( log(val) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int ii = 0; ii < w; ii++)
        {
            for(int jj = 0; jj < h; jj++)
            {
                val = 0;
                x = ii; y = (h-1) - jj;
                for(auto & grid : o)
                    val += (grid->point(ii+loc[0],jj+loc[1]) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);

    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (std::vector<cplx_grid_ptr> o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min, max;
    double temp1 = 0, temp2 = 0;
    for(auto & grid: o)
    {
        if (part == PLOTTYPE::REAL)
            std::tie(temp1, temp2) = findMinMaxReal(grid);
        else if (part == PLOTTYPE::IMAG)
            std::tie(temp1, temp2) = findMinMaxImag(grid);
        else if(part == PLOTTYPE::MAG || part == PLOTTYPE::POW)
            std::tie(temp1, temp2) = findMinMaxAbs(grid);
        min += temp1;
        max += temp2;
    }

    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    if(part == PLOTTYPE::POW)
    {
        max *= max;
        min *= min;
    }
    double diff = max - min;

    // int r,g,b;
    if (part == PLOTTYPE::REAL)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0.0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += (grid->point(i+loc[0],j+loc[1]).real() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::IMAG)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0.0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += (grid->point(i+loc[0],j+loc[1]).imag() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::MAG)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0.0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += sqrt( std::real( grid->point(i+loc[0],j+loc[1]) * std::conj(grid->point(i+loc[0],j+loc[1]) ) ) );
                val = (val - min) / diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::POW)
    {
        int x = 0, y = 0;
        double val = 0.0;
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                val = 0.0;
                x = i; y = (h-1) - j;
                for(auto& grid : o)
                    val += std::real( grid->point(i+loc[0],j+loc[1]) * std::conj(grid->point(i+loc[0],j+loc[1]) ) );
                val = (val - min) / diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}

void GridToBitMap (int_grid_ptr o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);

    int min = 0, max = 0;
    std::tie(min,max) = findMinMaxReal(o);
    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    if(part == PLOTTYPE::POW)
    {
        min = 0.0;
        max = std::max(min*min, max*max);
    }
    double diff = max - min;
    if(part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (pow(o->point(i+loc[0],j+loc[1]), 2.0) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
        for(int ii = 0; ii < w; ii++)
        {
            for(int jj = 0; jj < h; jj++)
            {
                int x = ii, y = (h-1) - jj;
                double val = (o->point(ii+loc[0],jj+loc[1]) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);

    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (real_grid_ptr o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min = 0.0, max = 0.0;
    if(part == PLOTTYPE::POW)
    {
        std::tie(min, max) = findMinMaxPower(o);
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        std::tie(min, max) = findMinMaxPower(o);
        min = log(min);
        max = log(max);
    }
    else
    {
        std::tie(min, max) = findMinMaxReal(o);
    }

    double diff = max - min;
    if(part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (pow(o->point(i+loc[0],j+loc[1]), 2.0) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = ( log(pow(o->point(i+loc[0],j+loc[1]), 2.0) ) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i+loc[0],j+loc[1]) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (cplx_grid_ptr o, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min, max;

    if (part == PLOTTYPE::REAL)
        std::tie(min,max) = findMinMaxReal(o);
    else if (part == PLOTTYPE::IMAG)
        std::tie(min,max) = findMinMaxImag(o);
    else if(part == PLOTTYPE::MAG || part == PLOTTYPE::POW)
        std::tie(min,max) = findMinMaxAbs(o);

    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    if(part == PLOTTYPE::POW)
    {
        max *= max;
        min *= min;
    }
    double diff = max - min;
    if (diff == 0.0)
    {
        max  = 1.0;
        min  = 0.0;
        diff = 1.0;
    }

    // int r,g,b;
    if (part == PLOTTYPE::REAL)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i+loc[0],j+loc[1]).real() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::IMAG)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i+loc[0],j+loc[1]).imag() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::MAG)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (sqrt(std::real(o->point(i+loc[0],j+loc[1]) * std::conj(o->point(i+loc[0],j+loc[1])))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (std::real(o->point(i+loc[0],j+loc[1]) * std::conj(o->point(i+loc[0],j+loc[1]))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (real_grid_ptr o_1, real_grid_ptr o_2, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min_1, min_2, max_1, max_2;
    std::tie(min_1,max_1) = findMinMaxReal(o_1);
    std::tie(min_2,max_2) = findMinMaxReal(o_2);
    double min = 0.0;
    double max = std::max(min_1*min_1 + min_2*min_2, max_1*max_1 + max_2*max_2);
    double diff = max - min;
    for(int i = 0; i < w; i++)
    {
        for(int j = 0; j < h; j++)
        {
            int x = i, y = (h-1) - j;
            double val = ((pow(o_1->point(i+loc[0],j+loc[1]), 2.0) + pow(o_2->point(i+loc[0],j+loc[1]), 2.0)) - min)/diff;
            img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
            img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
            img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (cplx_grid_ptr o_1, cplx_grid_ptr o_2, std::string filename, PLOTTYPE part, std::array<int,3> loc, std::array<int,3> sz)
{
    int w = sz[0];//o->x();
    int h = sz[1];//o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min_1, min_2, max_1, max_2;

    std::tie(min_1,max_1) = findMinMaxAbs(o_1);
    std::tie(min_2,max_2) = findMinMaxAbs(o_2);
    double min = 0.0;
    double max = std::max(min_1*min_1 + min_2*min_2, max_1*max_1 + max_2*max_2);

    double diff = max - min;
    if (diff == 0.0)
    {
        max  = 1.0;
        min  = 0.0;
        diff = 1.0;
    }

    // int r,g,b;
    else if (part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = ((std::real(o_1->point(i+loc[0],j+loc[1]) * std::conj(o_1->point(i+loc[0],j+loc[1]))) + std::real(o_2->point(i+loc[0],j+loc[1]) * std::conj(o_2->point(i+loc[0],j+loc[1])))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (int_grid_ptr o, std::string filename, PLOTTYPE part)
{
    int w = o->x();
    int h = o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);

    int min = 0, max = 0;
    std::tie(min,max) = findMinMaxReal(o);
    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    if(part == PLOTTYPE::POW)
    {
        min = 0.0;
        max = std::max(min*min, max*max);
    }
    double diff = max - min;
    if(part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (pow(o->point(i,j), 2.0) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i,j) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);

    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (real_grid_ptr o, std::string filename, PLOTTYPE part)
{
    int w = o->x();
    int h = o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min, max;
    std::tie(min,max) = findMinMaxReal(o);
    if(part == PLOTTYPE::POW)
    {
        std::tie(min, max) = findMinMaxPower(o);
        min = min*min;
        max = max*max;
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        std::tie(min, max) = findMinMaxPower(o);
        min = log(min);
        max = log(max);
    }
    else if(part == PLOTTYPE::MAG)
    {
        std::tie(min, max) = findMinMaxPower(o);
    }
    else
    {
        std::tie(min, max) = findMinMaxReal(o);
    }
    if(min == max)
    {
        min -= 0.05;
        max += 0.05;
    }
    double diff = max - min;

    if(part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (pow(o->point(i,j), 2.0) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if(part == PLOTTYPE::LNPOW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = ( log(pow(o->point(i,j), 2.0) ) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i,j) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (cplx_grid_ptr o, std::string filename, PLOTTYPE part)
{
    int w = o->x();
    int h = o->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min, max;

    if (part == PLOTTYPE::REAL)
        std::tie(min,max) = findMinMaxReal(o);
    else if (part == PLOTTYPE::IMAG)
        std::tie(min,max) = findMinMaxImag(o);
    else if(part == PLOTTYPE::MAG || part == PLOTTYPE::POW)
        std::tie(min,max) = findMinMaxAbs(o);

    if (min == max)
    {
        max *= 1.1;
        min *=  0.9;
    }
    if(part == PLOTTYPE::POW)
    {
        max *= max;
        min *= min;
    }
    double diff = max - min;
    if (diff == 0.0)
    {
        max  = 1.0;
        min  = 0.0;
        diff = 1.0;
    }

    // int r,g,b;
    if (part == PLOTTYPE::REAL)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i,j).real() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::IMAG)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (o->point(i,j).imag() - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::MAG)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (sqrt(std::real(o->point(i,j) * std::conj(o->point(i,j)))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }
    else if (part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = (std::real(o->point(i,j) * std::conj(o->point(i,j))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (real_grid_ptr o_1, real_grid_ptr o_2, std::string filename, PLOTTYPE part)
{
    int w = o_1->x();
    int h = o_1->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min_1, min_2, max_1, max_2;
    std::tie(min_1,max_1) = findMinMaxReal(o_1);
    std::tie(min_2,max_2) = findMinMaxReal(o_2);
    double min = 0.0;
    double max = std::max(min_1*min_1 + min_2*min_2, max_1*max_1 + max_2*max_2);
    double diff = max - min;
    for(int i = 0; i < w; i++)
    {
        for(int j = 0; j < h; j++)
        {
            int x = i, y = (h-1) - j;
            double val = ((pow(o_1->point(i,j), 2.0) + pow(o_2->point(i,j), 2.0)) - min)/diff;
            img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
            img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
            img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}
void GridToBitMap (cplx_grid_ptr o_1, cplx_grid_ptr o_2, std::string filename, PLOTTYPE part)
{
    int w = o_1->x();
    int h = o_1->y();
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

    std::vector<char> img(3*w*h);
    double min_1, min_2, max_1, max_2;

    std::tie(min_1,max_1) = findMinMaxAbs(o_1);
    std::tie(min_2,max_2) = findMinMaxAbs(o_2);
    double min = 0.0;
    double max = std::max(min_1*min_1 + min_2*min_2, max_1*max_1 + max_2*max_2);

    double diff = max - min;
    if (diff == 0.0)
    {
        max  = 1.0;
        min  = 0.0;
        diff = 1.0;
    }

    // int r,g,b;
    else if (part == PLOTTYPE::POW)
    {
        for(int i = 0; i < w; i++)
        {
            for(int j = 0; j < h; j++)
            {
                int x = i, y = (h-1) - j;
                double val = ((std::real(o_1->point(i,j) * std::conj(o_1->point(i,j))) + std::real(o_2->point(i,j) * std::conj(o_2->point(i,j)))) - min)/diff;
                img[(x+y*w)*3+2] = (unsigned char)(toRValue(val));
                img[(x+y*w)*3+1] = (unsigned char)(toGValue(val));
                img[(x+y*w)*3+0] = (unsigned char)(toBValue(val));
            }
        }
    }

    unsigned char bmpfileheader[14] = {'B','M',  0,0,0,0, 0,0,0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0,24,0};
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

    // boost::filesystem::path p(filename.c_str());
    // if( p.remove_filename().std::string() != "" )
    //     boost::filesystem::create_directories(p.remove_filename());
    std::ofstream outMap(filename, std::ios::out | std::ios::binary);
    outMap.write(reinterpret_cast<char *>(bmpfileheader),sizeof(bmpfileheader));
    outMap.write(reinterpret_cast<char *>(bmpinfoheader),sizeof(bmpinfoheader));
    for(int i=0; i<h; i++)
    {
        outMap.write(img.data()+(w*(h-i-1)*3),3*w);
        outMap.write(reinterpret_cast<char *>(bmppad),1*(4-(w*3)%4)%4);
    }
    outMap.close();
}


std::tuple<double,double> findMinMaxReal(real_grid_ptr &o)
{
    int maxLoc    = idamax_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
    double maxVal =  1.05 * std::abs( *(o->data() + maxLoc) );
    double minVal =  -1.0 * maxVal;
    return std::make_tuple(minVal,maxVal);
}

std::tuple<double,double> findMinMaxPower(real_grid_ptr &o)
{
//    int minLoc    = idamin_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
    int maxLoc    = idamax_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
    int minLoc    = idamin_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
    double maxVal =  1.01*std::abs(std::real( *(o->data() + maxLoc) ) );
    double minVal =  0.99*std::abs(std::real( *(o->data() + minLoc) ) ); // std::abs(std::real( *(o->data() + minLoc) ));
    return std::make_tuple(minVal,maxVal);
}

// std::tuple<double,double> findMinMaxPower(real_grid_ptr &o)
// {
//     int minLoc    = idamin_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
//     int maxLoc    = idamax_(o->size(),reinterpret_cast<double*>(o->data()), 1) - 1;
//     double maxVal =  1.01*std::abs(std::real( *(o->data() + maxLoc) ));
//     double minVal =  minLoc;
//     return std::make_tuple(minVal,maxVal);
// }

std::tuple<double,double> findMinMaxReal(cplx_grid_ptr &o)
{
//    int minLoc    = idamin_(o->size(),reinterpret_cast<double*>(o->data()), 2) - 1;
    int maxLoc    = idamax_(o->size(),reinterpret_cast<double*>(o->data()), 2) - 1;
    double maxVal =  1.05*std::abs(std::real( *(o->data() + maxLoc) ));
    double minVal =  -1.0 * maxVal;
    return std::make_tuple(minVal,maxVal);
}
std::tuple<int,int> findMinMaxReal(int_grid_ptr &o)
{
//    int minLoc    = isamin_(o->size(),reinterpret_cast<int*>(o->data()), 1) - 1;
    int maxLoc    = isamax_(o->size(),reinterpret_cast<int*>(o->data()), 1) - 1;
    int maxVal =  1.05*std::abs(std::real( *(o->data() + maxLoc) ));
    int minVal =  -1.0 * maxVal;
    return std::make_tuple(minVal,maxVal);
}
std::tuple<double,double> findMinMaxImag(cplx_grid_ptr &o)
{
//    int minLoc    = idamin_(o->size(),reinterpret_cast<double*>(o->data()) + 1, 2) - 1;
    int maxLoc    = idamax_(o->size(),reinterpret_cast<double*>(o->data()) + 1, 2) - 1;
    double maxVal = 1.05*std::abs(imag( *(o->data() + maxLoc) ));
    double minVal = -1.0 * maxVal;
    return std::make_tuple(minVal,maxVal);
}

std::tuple<double,double> findMinMaxAbs(cplx_grid_ptr &o)
{
//  int minLoc    = izamin_(o->size(), o->data(), 1) - 1;
  int maxLoc    = izamax_(o->size(), o->data(), 1) - 1;
  double maxVal = 1.05*std::abs( *(o->data() + maxLoc) );
  double minVal = -1.0 * maxVal;
  return std::make_tuple(minVal,maxVal);
}

