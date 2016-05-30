#pragma once
#ifndef _FRAME_INCLUDE_ONCE
#define _FRAME_INCLUDE_ONCE

#include "math.hpp"
#include <vector>
#include <fstream>

class Frame
{
public:

	// Constructor. Reserves the memory for the frame.
	Frame(const Vec2i& screenResolution) : _ScreenResolution(screenResolution)
	{
		_AccumColor.resize(screenResolution.x*screenResolution.y, Vec3d(0,0,0));
		_AccumCount.resize(screenResolution.x*screenResolution.y, 0);
	}

	// Adds a pixel measurement to the statistics.
	void AddPixel(const Vec2i& pxCoord, const Vec3d& color)
	{
		int linearIndex = pxCoord.y * _ScreenResolution.x + pxCoord.x;
		_AccumColor[linearIndex] += color;
		_AccumCount[linearIndex] += 1;
	}

	// Exports the current frame to a pfm file.
	void ExportPfm(const char *filename)
	{
		// create and replace the file
		std::ofstream headerWriter(filename);
		if (headerWriter.is_open())
		{
			// write the header
			headerWriter << "PF" << std::endl;
			headerWriter << _ScreenResolution.x << " " << _ScreenResolution.y << std::endl;
			headerWriter.close();
		}

		// open the file in binary mode
		std::ofstream rawWriter(filename, std::ios::binary | std::ios::app);

		std::vector<float> data(_ScreenResolution.x * _ScreenResolution.y * 3);
		for (int i = 0; i < _ScreenResolution.x * _ScreenResolution.y; ++i) {
			Vec3d color = _AccumCount[i] == 0 ? Vec3d(1,1,1) : _AccumColor[i] / _AccumCount[i];
			data[3*i+0] = (float)color.x();
			data[3*i+1] = (float)color.y();
			data[3*i+2] = (float)color.z();
		}

		// write the image line by line (reverse order)
		int lSize = _ScreenResolution.x*sizeof(float)*3;
		//for (int y = _ScreenResolution.y - 1; y >= 0; y--)
		for (int y = 0; y < _ScreenResolution.y; ++y)
			rawWriter.write((char*)data.data() + lSize*y, lSize);
		rawWriter.close();
	}

	// Exports the current frame to a bmp file.
	void ExportBmp(const char* filename, double contrast, double brightness, double gamma)
	{
		unsigned char file[14] = {
			'B','M', // magic
			0,0,0,0, // size in bytes
			0,0, // app data
			0,0, // app data
			40+14,0,0,0 // start of data offset
		};
		unsigned char info[40] = {
			40,0,0,0, // info hd size
			0,0,0,0, // width
			0,0,0,0, // heigth
			1,0, // number color planes
			24,0, // bits per pixel
			0,0,0,0, // compression is none
			0,0,0,0, // image bits size
			0x13,0x0B,0,0, // horz resoluition in pixel / m
			0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
			0,0,0,0, // #colors in pallete
			0,0,0,0, // #important colors
			};

		int w=_ScreenResolution.x;
		int h=_ScreenResolution.y;

		int padSize  = (4-(w*3)%4)%4;
		int sizeData = w*h*3 + h*padSize;
		int sizeAll  = sizeData + sizeof(file) + sizeof(info);

		file[ 2] = (unsigned char)( sizeAll    );
		file[ 3] = (unsigned char)( sizeAll>> 8);
		file[ 4] = (unsigned char)( sizeAll>>16);
		file[ 5] = (unsigned char)( sizeAll>>24);

		info[ 4] = (unsigned char)( w   );
		info[ 5] = (unsigned char)( w>> 8);
		info[ 6] = (unsigned char)( w>>16);
		info[ 7] = (unsigned char)( w>>24);

		info[ 8] = (unsigned char)( h    );
		info[ 9] = (unsigned char)( h>> 8);
		info[10] = (unsigned char)( h>>16);
		info[11] = (unsigned char)( h>>24);

		info[20] = (unsigned char)( sizeData    );
		info[21] = (unsigned char)( sizeData>> 8);
		info[22] = (unsigned char)( sizeData>>16);
		info[23] = (unsigned char)( sizeData>>24);

		std::ofstream stream(filename, std::ios::binary | std::ios::out);
		stream.write( (char*)file, sizeof(file) );
		stream.write( (char*)info, sizeof(info) );

		unsigned char pad[3] = {0,0,0};

		for ( int y=0; y<h; y++ )
		{
			for ( int x=0; x<w; x++ )
			{
				Vec3d color = _AccumCount[y*w+x] == 0 ? Vec3d(1,1,1) : _AccumColor[y*w+x] / _AccumCount[y*w+x];
				color = pow(max(0.f,(color-0.5)*contrast+0.5+brightness), 1.0/gamma);
				unsigned char pixel[3];
				pixel[0] = (unsigned char)(std::min(std::max(0.f, (float)color.z()), 1.f)*255);
				pixel[1] = (unsigned char)(std::min(std::max(0.f, (float)color.y()), 1.f)*255);
				pixel[2] = (unsigned char)(std::min(std::max(0.f, (float)color.x()), 1.f)*255);

				stream.write( (char*)pixel, 3 );
			}
			stream.write( (char*)pad, padSize );
		}
	}
	
private:

	std::vector<Vec3d> _AccumColor;
	std::vector<double> _AccumCount;
	Vec2i _ScreenResolution;			// resolution of the viewport
};

#endif