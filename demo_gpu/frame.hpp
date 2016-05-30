#pragma once
#ifndef _FRAME_INCLUDE_ONCE
#define _FRAME_INCLUDE_ONCE

#include "math.hpp"
#include <vector>
#include <random>
#include <fstream>
#include <d3d11.h>

class Frame
{
public:

	// Constructor. Creates the resources used to store the pixel statistics on both the CPU and GPU.
	Frame(const Vec2i& screenResolution, ID3D11Device* device) : 
		_ScreenResolution(screenResolution), 
		_BufAccumRadiance(NULL), 
		_BufAccumCount(NULL),
		_BufAccumRadianceStaging(NULL),
		_BufAccumCountStaging(NULL),
		_UavAccumRadiance(NULL),
		_UavAccumCount(NULL),
		_BufRandomSeeds(NULL),
		_UavRandomSeeds(NULL)
	{
		_AccumRadiance.resize(screenResolution.x*screenResolution.y, Vec3f(0,0,0));
		_AccumCount.resize(screenResolution.x*screenResolution.y, 0);

		{
			D3D11_SUBRESOURCE_DATA init;
			ZeroMemory(&init, sizeof(D3D11_SUBRESOURCE_DATA));
			init.pSysMem = _AccumRadiance.data();

			D3D11_BUFFER_DESC bufDesc;
			ZeroMemory(&bufDesc, sizeof(D3D11_BUFFER_DESC));
			bufDesc.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
			bufDesc.ByteWidth = screenResolution.x * screenResolution.y * sizeof(float)*3;
			bufDesc.StructureByteStride = sizeof(float)*3;
			bufDesc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
			bufDesc.Usage = D3D11_USAGE_DEFAULT;
			bufDesc.CPUAccessFlags = 0;
			device->CreateBuffer(&bufDesc, &init, &_BufAccumRadiance);
			device->CreateUnorderedAccessView(_BufAccumRadiance, NULL, &_UavAccumRadiance);

			bufDesc.BindFlags = 0;
			bufDesc.Usage = D3D11_USAGE_STAGING;
			bufDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
			device->CreateBuffer(&bufDesc, NULL, &_BufAccumRadianceStaging);
		}

		{
			D3D11_SUBRESOURCE_DATA init;
			ZeroMemory(&init, sizeof(D3D11_SUBRESOURCE_DATA));
			init.pSysMem = &_AccumCount[0];

			D3D11_BUFFER_DESC bufDesc;
			ZeroMemory(&bufDesc, sizeof(D3D11_BUFFER_DESC));
			bufDesc.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
			bufDesc.ByteWidth = screenResolution.x * screenResolution.y * sizeof(float);
			bufDesc.StructureByteStride = sizeof(float);
			bufDesc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
			bufDesc.Usage = D3D11_USAGE_DEFAULT;
			bufDesc.CPUAccessFlags = 0;
			device->CreateBuffer(&bufDesc, &init, &_BufAccumCount);
			device->CreateUnorderedAccessView(_BufAccumCount, NULL, &_UavAccumCount);

			bufDesc.BindFlags = 0;
			bufDesc.Usage = D3D11_USAGE_STAGING;
			bufDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
			device->CreateBuffer(&bufDesc, NULL, &_BufAccumCountStaging);
		}

		{
			D3D11_BUFFER_DESC bufDesc;
			ZeroMemory(&bufDesc, sizeof(D3D11_BUFFER_DESC));
			bufDesc.BindFlags = D3D11_BIND_UNORDERED_ACCESS;
			bufDesc.ByteWidth = screenResolution.x * screenResolution.y * sizeof(int)*2;
			bufDesc.StructureByteStride = sizeof(int)*2;
			bufDesc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;

			std::default_random_engine rng;
			std::uniform_int_distribution<int> rndi;
		
			std::vector<int> seeds(screenResolution.x * screenResolution.y*2);
			for (size_t i=0; i<seeds.size(); ++i)
				seeds[i] = rndi(rng);

			D3D11_SUBRESOURCE_DATA init;
			ZeroMemory(&init, sizeof(D3D11_SUBRESOURCE_DATA));
			init.pSysMem = &seeds[0];

			device->CreateBuffer(&bufDesc, &init, &_BufRandomSeeds);
			device->CreateUnorderedAccessView(_BufRandomSeeds, NULL, &_UavRandomSeeds);
		}
	}

	~Frame()
	{
		if (_BufAccumRadiance) _BufAccumRadiance->Release();
		if (_BufAccumCount) _BufAccumCount->Release();
		if (_BufRandomSeeds) _BufRandomSeeds->Release();
		if (_BufAccumRadianceStaging) _BufAccumRadianceStaging->Release();
		if (_BufAccumCountStaging) _BufAccumCountStaging->Release();
		if (_UavAccumRadiance) _UavAccumRadiance->Release();
		if (_UavAccumCount) _UavAccumCount->Release();
		if (_UavRandomSeeds) _UavRandomSeeds->Release();
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
			Vec3f color = _AccumCount[i] == 0 ? Vec3f(1,1,1) : _AccumRadiance[i] / _AccumCount[i];
			data[3*i+0] = color.x;
			data[3*i+1] = color.y;
			data[3*i+2] = color.z;
		}

		// write the image line by line (reverse order)
		int lSize = _ScreenResolution.x*sizeof(float)*3;
		//for (int y = _ScreenResolution.y - 1; y >= 0; y--)
		for (int y = 0; y < _ScreenResolution.y; ++y)
			rawWriter.write((char*)data.data() + lSize*y, lSize);
		rawWriter.close();
	}

	// Exports the current frame to a bmp file.
	void ExportBmp(const char* filename, float contrast, float brightness, float gamma)
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
				Vec3f color = _AccumCount[y*w+x] == 0 ? Vec3f(1,1,1) : _AccumRadiance[y*w+x] / _AccumCount[y*w+x];
				color = pow(max(0.f,(color-0.5f)*contrast+0.5f+brightness), 1.0f/gamma);
				unsigned char pixel[3];
				pixel[0] = (unsigned char)(std::min(std::max(0.f, color.z), 1.f)*255);
				pixel[1] = (unsigned char)(std::min(std::max(0.f, color.y), 1.f)*255);
				pixel[2] = (unsigned char)(std::min(std::max(0.f, color.x), 1.f)*255);

				stream.write( (char*)pixel, 3 );
			}
			stream.write( (char*)pad, padSize );
		}
	}
	
	ID3D11UnorderedAccessView* GetUavAccumRadiance() { return _UavAccumRadiance; }
	ID3D11UnorderedAccessView* GetUavAccumCount() { return _UavAccumCount; }
	ID3D11UnorderedAccessView* GetUavRandomSeeds() { return _UavRandomSeeds; }

	// Copies the pixel statistics from the GPU to the CPU via staging buffers
	void CopyToHost(ID3D11DeviceContext* immediateContext)
	{
		immediateContext->CopyResource(_BufAccumRadianceStaging, _BufAccumRadiance);
		immediateContext->CopyResource(_BufAccumCountStaging, _BufAccumCount);

		D3D11_MAPPED_SUBRESOURCE mapped;
		if (SUCCEEDED(immediateContext->Map(_BufAccumRadianceStaging, 0, D3D11_MAP_READ, 0, &mapped)))
		{
			int sizeinbyte = _ScreenResolution.x*_ScreenResolution.y*sizeof(Vec3f);
			memcpy_s(_AccumRadiance.data(), sizeinbyte, mapped.pData, sizeinbyte);
			immediateContext->Unmap(_BufAccumRadianceStaging, 0);
		}

		if (SUCCEEDED(immediateContext->Map(_BufAccumCountStaging, 0, D3D11_MAP_READ, 0, &mapped)))
		{
			int sizeinbyte = _ScreenResolution.x*_ScreenResolution.y*sizeof(float);
			memcpy_s(_AccumCount.data(), sizeinbyte, mapped.pData, sizeinbyte);
			immediateContext->Unmap(_BufAccumCountStaging, 0);
		}
	}

private:

	Vec2i _ScreenResolution;				// resolution of the viewport

	std::vector<Vec3f> _AccumRadiance;		// accumulates the radiance (CPU-side)
	std::vector<float> _AccumCount;			// counter for the number of measurements per pixel (CPU-side)
	
	ID3D11Buffer* _BufAccumRadiance;		// accumulates the radiance (GPU-side), float3 structured buffer
	ID3D11Buffer* _BufAccumCount;			// counter for the number of measurements per pixel (GPU-side), float structured buffer

	ID3D11Buffer* _BufAccumRadianceStaging;	// staging buffer that is used to copy _BufAccumRadiance (GPU) to _AccumRadiance (CPU)
	ID3D11Buffer* _BufAccumCountStaging;	// staging buffer that is used to copy _BufAccumRadiance (GPU) to _AccumRadiance (CPU)

	ID3D11UnorderedAccessView* _UavAccumRadiance;	// Unordered access view that gives compute shaders read/write access to _BufAccumRadiance.
	ID3D11UnorderedAccessView* _UavAccumCount;		// Unordered access view that gives compute shaders read/write access to _BufAccumCount.

	ID3D11Buffer* _BufRandomSeeds;					// random seeds per pixel (GPU), int2
	ID3D11UnorderedAccessView* _UavRandomSeeds;		// Unordered access view that gives compute shaders read/write access to _BufRandomSeeds.

};

#endif