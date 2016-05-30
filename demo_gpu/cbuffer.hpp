#pragma once

#include <D3D11.h>

// This class creates a Direct3D 11 constant buffer

template <typename T>
class ConstantBuffer
{
	public:

		// Constructor.
		ConstantBuffer() : Data(), _Buffer(NULL) {}
		// Destructor. Releases the Direct3D resource.
		~ConstantBuffer() { Release(); }

		// Copy of the data on the CPU side
		T Data;

		// Gets the Direct3D resource.
		ID3D11Buffer* GetBuffer() const { return _Buffer; }

		// Creates the Direct3D resource
		bool Create(ID3D11Device* Device)
		{
			int size = sizeof(T);
			if ((size & 15) != 0)	// take care of padding
			{
				size >>= 4;
				size++;
				size <<= 4;
			}

			D3D11_BUFFER_DESC desc;	
			ZeroMemory(&desc, sizeof(D3D11_BUFFER_DESC));
			desc.Usage = D3D11_USAGE_DYNAMIC;
			desc.ByteWidth = size;
			desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
			desc.MiscFlags = 0;

			D3D11_SUBRESOURCE_DATA init;
			ZeroMemory(&init, sizeof(D3D11_SUBRESOURCE_DATA));
			init.pSysMem = &Data;
			if (S_OK != Device->CreateBuffer(&desc, &init, &_Buffer))
				return false;

			return true;
		}

		// Releases the Direct3D resource
		void Release()
		{
			if (_Buffer) { _Buffer->Release(); _Buffer = NULL; }
		}

		// Maps the CPU content to the GPU
		void UpdateBuffer(ID3D11DeviceContext* ImmediateContext)
		{
			if (!_Buffer) return;
			D3D11_MAPPED_SUBRESOURCE MappedSubResource;		
			ImmediateContext->Map(_Buffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedSubResource );			
			*(T*)MappedSubResource.pData = Data;
			ImmediateContext->Unmap(_Buffer, 0 );		
		}

	private:

		// Pointer to the Direct3D Resource
		ID3D11Buffer* _Buffer;
		
};