#ifndef _D3D_INCLUDE_ONCE
#define _D3D_INCLUDE_ONCE

#include <dxgi.h>
#include <d3d11_1.h>
#include <fstream>

class D3D
{
public:

	// Constructor. Creates a Direct3D device.
	D3D() : _Device(NULL), _ImmediateContext(NULL)
	{
		IDXGIFactory* factory;
		if (FAILED(CreateDXGIFactory(__uuidof(IDXGIFactory), (void**)&factory))) return;
		
		HRESULT hr = S_OK;

		UINT createDeviceFlags = 0;
#ifdef _DEBUG
		createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif

		D3D_DRIVER_TYPE driverTypes[] =
		{
			D3D_DRIVER_TYPE_UNKNOWN,
		};
		UINT numDriverTypes = sizeof(driverTypes) / sizeof(driverTypes[0]);

		D3D_FEATURE_LEVEL featureLevels[] = { D3D_FEATURE_LEVEL_11_0 };
		unsigned int numFeatureLevels = 1;
		D3D_FEATURE_LEVEL usedFeatureLevel = D3D_FEATURE_LEVEL_11_0;

		// iterate the display adapters and look for a DirectX 11 capable device.
		for (UINT driverTypeIndex = 0; driverTypeIndex < numDriverTypes; driverTypeIndex++)
		{
			D3D_DRIVER_TYPE driverType = driverTypes[driverTypeIndex];

			IDXGIAdapter* adapter;
			for (UINT i = 0;
				factory->EnumAdapters(i, &adapter) != DXGI_ERROR_NOT_FOUND;
				++i)
			{
				DXGI_ADAPTER_DESC desc;
				adapter->GetDesc(&desc);

				std::wstring adapterDesc(desc.Description);
				hr = D3D11CreateDevice(adapter, driverType, NULL, createDeviceFlags, featureLevels, numFeatureLevels, D3D11_SDK_VERSION, &_Device, &usedFeatureLevel, &_ImmediateContext);
				if (SUCCEEDED(hr))
				{
					wprintf(L"D3D is using: %s\n", adapterDesc.c_str());
					break;
				}
			}

			if (!_Device) {
				printf("Couldn't create DirectX device.\nMaybe DirectX 11 is not supported on your machine?\n");
				exit(-1);
			}

			if (adapter) adapter->Release();

			if (SUCCEEDED(hr))
				break;
		}
		factory->Release();
	}

	// Destructor. Releases the resources.
	~D3D()
	{
		_Device->Release();
		_ImmediateContext->Release();
	}

	// Utility function that loads a compute shader from file
	bool LoadComputeShaderFromFile(const char* path, ID3D11ComputeShader** outShader)
	{
		unsigned int size;
		char* memblock = NULL;
		std::ifstream file(path, std::ios::in | std::ios::binary | std::ios::ate);
		if (file.is_open())
		{
			size = (unsigned int)file.tellg();
			memblock = new char[size];
			file.seekg(0, std::ios::beg);
			file.read(memblock, size);
			file.close();

			if (FAILED(_Device->CreateComputeShader(memblock, size, NULL, outShader))) {
				delete[] memblock;
				return false;
			}
			delete[] memblock;
		}
		return true;
	}

	// Gets the Direct3D device.
	ID3D11Device* GetDevice() { return _Device; }
	// Gets the immediate context.
	ID3D11DeviceContext* GetImmediateContext() { return _ImmediateContext; }

private:

	ID3D11Device* _Device;						// Direct3D device
	ID3D11DeviceContext* _ImmediateContext;		// Immediate context
};

#endif