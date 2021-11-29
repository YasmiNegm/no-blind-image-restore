#include <iostream>
#include <complex>
#include <fstream>
#include <cstring>
#include <windows.h>

#define pi 3.1415926536

using namespace std;

#pragma pack(1)

typedef struct RGB
{
	double R, G, B;

};

typedef struct cRGB
{
	complex<double> R, G, B;

};

struct f_info
{
	unsigned char signature[2];	
	unsigned int sizefile;		
	unsigned int reserved;		
	unsigned int addr_offset;	
};

struct pic_info
{
	unsigned int Size;        	
	unsigned int M;       	
	unsigned int N;      	
	unsigned short int  Planes;   
	unsigned short int  BitCount; 
	unsigned int Compression; 	
	unsigned int SizeImage;   	
	unsigned int XPelsPerMeter;	
	unsigned int YPelsPerMeter;	
	unsigned int ClrUsed;     	
	unsigned int ClrImportant;	
};

void filtrBFHF(cRGB* Iimage, cRGB* (&Iimage1), int H, int W)
{
	int m = 0, h = 1, w = 1;
	for (int tmp_size = H; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, h *= 2);
	if (h < H) h *= 2;
	m = 0;
	for (int tmp_size = W; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, w *= 2);
	if (w < W) w *= 2;
	if (w < h) w = h;
	Iimage1 = new cRGB[w * w];
	int CY = w / 2;
	int CX = w / 2;
	int n = 2;
	double Dp = sqrt(CX * CX + CY * CY) * 2;
	for (int u = 0; u < w; u++) {		//H
		for (int v = 0; v < w; v++) {	//W
			double D = sqrt((u - CX) * (u - CX) + (v - CY) * (v - CY));
			double Huv = 1 / (1 + pow(Dp / D, 2 * n));
			Iimage1[u * w + v].R = Iimage[u * w + v].R * Huv;
			Iimage1[u * w + v].G = Iimage[u * w + v].G * Huv;
			Iimage1[u * w + v].B = Iimage[u * w + v].B * Huv;
		}
	}
}


void filtrWnr(cRGB* Iimage, cRGB* (&Iimage1), int H, int W)
{
	int m = 0, h = 1, w = 1;
	for (int tmp_size = H; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, h *= 2);
	if (h < H) h *= 2;
	m = 0;

	for (int tmp_size = W; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, w *= 2);
	if (w < W) w *= 2;
	if (w < h) w = h;

	Iimage1 = new cRGB[w * w];
	int CY = w / 2;
	int CX = w / 2;

	int n = 1;
	double Dp = sqrt(CX * CX + CY * CY) * 0.08;
	for (int u = 0; u < w; u++) {		
		for (int v = 0; v < w; v++) {	
			double D = sqrt((u - CX) * (u - CX) + (v - CY) * (v - CY));
			double Huv = 1 / (1 + pow(D / Dp, 2 * n)); 

			Iimage1[u * w + v].R = (1 / Huv * (pow(abs(Huv),2) / (pow(abs(Huv), 2) + 0.01) ) ) *Iimage[u * w + v].R;
			Iimage1[u * w + v].G = (1 / Huv * (pow(abs(Huv), 2) / (pow(abs(Huv), 2) + 0.01) ) ) *Iimage[u * w + v].G;
			Iimage1[u * w + v].B = (1 / Huv * (pow(abs(Huv), 2) / (pow(abs(Huv), 2) + 0.01) ) ) *Iimage[u * w + v].B;
		}
	}
}

void BitReversing(cRGB* inputSignal, cRGB* outputSignal, int size, int step, int offset)
{
	int m = 0, n;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);
	n = 1 << m;
	for (unsigned int i = 0; i < n; ++i) {
		int k = 0, x = i;
		for (int j = 0; j < m; j++) {
			k <<= 1;
			k |= (x & 1);
			x >>= 1;
		}
		outputSignal[k * step + offset] = inputSignal[i * step + offset];
	}
}

void fft(cRGB* b, int size, int step, int offset) {
	int log2n = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, log2n++);
	const std::complex<double> J(0, 1);
	int n = 1 << log2n;
	for (int s = 1; s <= log2n; ++s) {
		int m = 1 << s;
		int m2 = m >> 1;
		std::complex<double> w(1, 0);
		std::complex<double> wm = exp(-J * (pi / m2));
		for (int j = 0; j < m2; ++j) {
			for (int k = j; k < n; k += m) {
				int K = k * step + offset;
				int K1 = (k + m2) * step + offset;
				std::complex<double> tR = w * b[K1].R;
				std::complex<double> tG = w * b[K1].G;
				std::complex<double> tB = w * b[K1].B;
				std::complex<double> uR = b[K].R;
				std::complex<double> uG = b[K].G;
				std::complex<double> uB = b[K].B;
				b[K].R = (uR + tR) / 2.0;
				b[K].G = (uG + tG) / 2.0;
				b[K].B = (uB + tB) / 2.0;
				b[K1].R = (uR - tR) / 2.0;
				b[K1].G = (uG - tG) / 2.0;
				b[K1].B = (uB - tB) / 2.0;
			}
			w *= wm;
		}
	}
}

void fft1d(cRGB* a, cRGB* b, int size, int step, int offset)
{
	BitReversing(a, b, size, step, offset);
	fft(b, size, step, offset);
}

void ttf2D(RGB* masin, cRGB* (&masout), int h1, int w1)//получение комплексного представления RGB
{
	int m = 0, h = 1, w = 1, i, j;
	for (int tmp_size = h1; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, h *= 2);
	if (h < h1) h *= 2;
	m = 0;
	for (int tmp_size = w1; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, w *= 2);
	if (w < w1) w *= 2;
	if (w < h) w = h;
	else
		if (w > h)  h = w;
	cRGB* cmasin = new cRGB[h * w];
	cRGB* ctemp = new cRGB[h * w];
	masout = new cRGB[h * w];
	double k = -1.0;
	for (i = 0; i < h1; i++) {
		k *= -1.0;
		double kk = k;
		for (j = 0; j < w1; j++) {
			cmasin[i * w + j].R.real(masin[i * w1 + j].R * kk);
			cmasin[i * w + j].R.imag(0);
			cmasin[i * w + j].G.real(masin[i * w1 + j].G * kk);
			cmasin[i * w + j].G.imag(0);
			cmasin[i * w + j].B.real(masin[i * w1 + j].B * kk);
			cmasin[i * w + j].B.imag(0);
			kk *= -1;
		}
		for (; j < w; j++) {
			cmasin[i * w + j].R.real(0);
			cmasin[i * w + j].R.imag(0);
			cmasin[i * w + j].G.real(0);
			cmasin[i * w + j].G.imag(0);
			cmasin[i * w + j].B.real(0);
			cmasin[i * w + j].B.imag(0);
		}
	}
	for (; i < h; i++)
		for (j = 0; j < w; j++) {
			cmasin[i * w + j].R.real(0);
			cmasin[i * w + j].R.imag(0);
			cmasin[i * w + j].G.real(0);
			cmasin[i * w + j].G.imag(0);
			cmasin[i * w + j].B.real(0);
			cmasin[i * w + j].B.imag(0);
		}
	for (int i = 0; i < w; i++)
		fft1d(cmasin, ctemp, h, w, i);
	for (int j = 0; j < h; j++)
		fft1d(ctemp, masout, w, 1, h * j);
	delete[] cmasin;
	delete[] ctemp;
}

void ifft(cRGB* b, int size, int step, int offset)
{
	const std::complex<double> J(0, 1);
	int log2n = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, log2n++);
	int n = 1 << log2n;
	for (int s = 1; s <= log2n; ++s) {
		int m = 1 << s;
		int m2 = m >> 1;
		std::complex<double> w(1, 0);
		std::complex<double> wm = exp(J * (pi / m2));
		for (int j = 0; j < m2; ++j) {
			for (int k = j; k < n; k += m) {
				std::complex<double> tR = w * b[(k + m2) * step + offset].R;
				std::complex<double> tG = w * b[(k + m2) * step + offset].G;
				std::complex<double> tB = w * b[(k + m2) * step + offset].B;
				std::complex<double> uR = b[k * step + offset].R;
				std::complex<double> uG = b[k * step + offset].G;
				std::complex<double> uB = b[k * step + offset].B;
				b[k * step + offset].R = (uR + tR);
				b[k * step + offset].G = (uG + tG);
				b[k * step + offset].B = (uB + tB);
				b[(k + m2) * step + offset].R = (uR - tR);
				b[(k + m2) * step + offset].G = (uG - tG);
				b[(k + m2) * step + offset].B = (uB - tB);
			}
			w *= wm;
		}
	}
}

void ifft1d(cRGB* a, cRGB* b, int size, int step, int offset)
{
	BitReversing(a, b, size, step, offset);
	ifft(b, size, step, offset);
}

void ittf2D(cRGB* cmasin, RGB* masout, int h1, int w1) 
{
	int m = 0, h = 1, w = 1;// , i, j;
	for (int tmp_size = h1; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, h *= 2);
	if (h < h1) h *= 2;
	m = 0;
	for (int tmp_size = w1; tmp_size > 1; tmp_size /= 2, m++);
	for (int i = 0; i < m; i++, w *= 2);
	if (w < w1) w *= 2;
	if (w < h) w = h;
	else
		if (w > h)  h = w;
	cRGB* cmasout = new cRGB[h * w];
	cRGB* ctemp = new cRGB[h * w];
	for (int i = 0; i < w; i++)
		ifft1d(cmasin, ctemp, h, w, i);
	for (int j = 0; j < h; j++)
		ifft1d(ctemp, cmasout, w, 1, h * j);
	double k = -1;
	for (int i = 0; i < h1; i++) {
		k *= -1;
		double kk = k;
		for (int j = 0; j < w1; j++) {
			masout[i * w1 + j].R = std::abs(cmasout[i * w + j].R);
			masout[i * w1 + j].G = std::abs(cmasout[i * w + j].G);
			masout[i * w1 + j].B = std::abs(cmasout[i * w + j].B);
			kk *= -1;
		}
	}
	delete[] cmasout;
	delete[] ctemp;
}


int main(int argc, char* argv[])
{
	char infile[100] = "deform.bmp";
	char outfile[100] = "repair.bmp";
	

	ifstream in(infile, std::ios::in | std::ios::binary);
	if (!in) {
		cout << "File not found\n";
		exit(-1);
	}


	struct f_info f_i;
	struct pic_info pic_i;

	in.read((char*)&f_i, sizeof(f_info));
	in.read((char*)&pic_i, sizeof(pic_info));

	int ln_str = (f_i.sizefile - 54) / pic_i.N; 
	RGB* image = new RGB[pic_i.N * pic_i.M]; 

	std::cout << "Picture size:\n";
	std::cout << "Width - " << pic_i.M << ", Height - " << pic_i.N << std::endl;

	for (auto i = 0; i < pic_i.N; i++)
	{
		for (auto j = 0; j < pic_i.M; j++) {

			char c;
			in.get(c);
			image[i * pic_i.M + j].B = (unsigned char)c;
			in.get(c);
			image[i * pic_i.M + j].G = (unsigned char)c;
			in.get(c);
			image[i * pic_i.M + j].R = (unsigned char)c;

		}

		for (auto k = 0; k < (ln_str - pic_i.M * 3); k++)
		{
			char d;
			in.get(d);
		}
	}
	in.close();

	
	//восстановление изображения

	cRGB* Iimage = nullptr;
	cRGB* Iimage2 = nullptr;
	cRGB* Iimage3 = nullptr;
	RGB* Iimage1 = new RGB[pic_i.N * pic_i.M];
	

	LARGE_INTEGER time1, time2;
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq); 
	QueryPerformanceCounter(&time1); 


	ttf2D(image, Iimage, pic_i.N, pic_i.M); 
	filtrWnr(Iimage, Iimage2, pic_i.N, pic_i.M);
	filtrWnr(Iimage2, Iimage3, pic_i.N, pic_i.M);
    ittf2D(Iimage3, Iimage1, pic_i.N, pic_i.M); 


	QueryPerformanceCounter(&time2);
	time2.QuadPart -= time1.QuadPart;
	double span = (double)time2.QuadPart / freq.QuadPart;
	std::cout << "Elasped time: " << span << std::endl;

	
	ofstream out(outfile, std::ios::out | std::ios::binary);
	out.write((char*)&f_i, sizeof(f_info));
	out.write((char*)&pic_i, sizeof(pic_info));

	for (auto i = 0; i < pic_i.N; i++) {
		for (auto j = 0; j < pic_i.M; j++) {
			out.put((char)(image[i * pic_i.M + j].B - Iimage1[i * pic_i.M + j].B ));
			out.put((char)(image[i * pic_i.M + j].G - Iimage1[i * pic_i.M + j].G ));
			out.put((char)(image[i * pic_i.M + j].R - Iimage1[i * pic_i.M + j].R ));
		}
		for (auto k = 0; k < (ln_str - pic_i.M * 3); k++) out.put(0);
	}

	out.close();
	delete[] image;
	delete[] Iimage;
	delete[] Iimage1;
}
