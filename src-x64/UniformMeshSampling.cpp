#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>

using namespace arma;

//   -*-   -*-   -*-

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

//   -*-   -*-   -*-

class CVec3Df {
    public:
        float X, Y, Z;
        
        CVec3Df() : X(0.0f), Y(0.0f), Z(0.0f) {
        }
        
        CVec3Df(float X, float Y, float Z) : X(X), Y(Y), Z(Z) {
        }
};

class CFace {
    public:
        int P1;
        int P2;
        int P3;
};

//   -*-   -*-   -*-

unsigned seedt = 0;

float RandomFloatt() {
    float val;
    *((unsigned *)&val) = (seedt & 8388607) | 1065353216;
    seedt = (seedt * 1664525) + 1013904223;
    if (val - 1.0f == 0.0f) {
        *((unsigned *)&val) = (seedt & 8388607) | 1065353216;
        seedt = (seedt * 1664525) + 1013904223;
    }
    return val - 1.0f;
}

CVec3Df SampleTriUniform(std::vector<CVec3Df> &vertices, CFace &f) {
	float u = RandomFloatt();
	float v = RandomFloatt();
	while (u + v > 1.0f) v = RandomFloatt();
	float w = 1.0f - u - v;
	CVec3Df P1 = vertices[f.P1];
	CVec3Df P2 = vertices[f.P2];
	CVec3Df P3 = vertices[f.P3];
	return CVec3Df(
	    (P1.X * u) + (P2.X * v) + (P3.X * w),
	    (P1.Y * u) + (P2.Y * v) + (P3.Y * w),
	    (P1.Z * u) + (P2.Z * v) + (P3.Z * w)
	);
}

int InverseTransformSampling(float *cdf, int size) {
    float val = RandomFloatt();
    int L = 0;
    int U = size - 1;
    while (L < U) {
        int M = (L + U) >> 1;
        if (cdf[M] < val) L = M + 1;
        else
            U = M; 
    }
    return L;
}

//   -*-   -*-   -*-

class CStringTokenizer {
	private:
		bool occArr[256];
		char lastDel;
		char *str;
		int ind;
		int len;
	public:
		CStringTokenizer(char *str) : str(str), ind(0), len(strlen(str)) {
			memset(occArr, 0, sizeof(bool) * 256);
		}

		char *NextTokenIE(char *delims);
		char *NextTokenAE(char *delims);
		void Dispose();
};

char *CStringTokenizer::NextTokenIE(char *delims) {
	char *res;
	int delLen = strlen(delims);
	
	if (ind > 0) str[ind - 1] = lastDel;
	if (ind < len) {
		for (int i = 0; i < delLen; ++i) occArr[delims[i]] = true;
		while ((ind < len) && (occArr[str[ind]])) ++ind;
		if (ind < len) {
			res = &str[ind];
			while ((ind < len) && (!occArr[str[ind]])) ++ind;
			lastDel = str[ind];
			str[ind++] = 0;
		} else
			res = NULL;
		for (int i = 0; i < delLen; ++i) occArr[delims[i]] = false;
	} else
		res = NULL;
	return res;
}

char *CStringTokenizer::NextTokenAE(char *delims) {
	char *res;
	int delLen = strlen(delims);
	
	if (ind > 0) str[ind - 1] = lastDel;
	if (ind < len) {
		res = &str[ind];
		for (int i = 0; i < delLen; ++i) occArr[delims[i]] = true;
		while ((ind < len) && (!occArr[str[ind]])) ++ind;
		lastDel = str[ind];
		str[ind++] = 0;
		for (int i = 0; i < delLen; ++i) occArr[delims[i]] = false;
	} else
		res = NULL;
	return res;
}

void CStringTokenizer::Dispose() {
	if (ind > 0) str[ind - 1] = lastDel;
}

//   -*-   -*-   -*-

// [[Rcpp::export]]
arma::mat SampleMeshUniform(std::string fName, int pointsNum) {
	char *buf;
	char *line;
	FILE *f = fopen(fName.c_str(), "rb+");
	float XMin = INFINITY, XMax = -INFINITY;
	float YMin = INFINITY, YMax = -INFINITY;
	float ZMin = INFINITY, ZMax = -INFINITY;
	int fSize;

	fseek(f, 0, SEEK_END);
	fSize = ftell(f);
	buf = (char *)malloc(sizeof(char) * (fSize + 1));
	fseek(f, 0, SEEK_SET);
	fread(buf, 1, fSize, f);
	buf[fSize] = 0;
	fclose(f);
	CStringTokenizer st1(buf);
	line = st1.NextTokenIE("\n\r");
	std::list<CVec3Df> vertices;
	std::list<CFace> faces;
	while (line != NULL) {
		char *token1;
		char *tokens1 = strdup(line);
		CStringTokenizer st2(tokens1);

		token1 = st2.NextTokenIE(" ");
		if (token1 != NULL) {
			if (strcmp(token1, "v") == 0) {
				CVec3Df v;

				sscanf(st2.NextTokenIE(" "), "%f", &v.X);
				sscanf(st2.NextTokenIE(" "), "%f", &v.Y);
				sscanf(st2.NextTokenIE(" "), "%f", &v.Z);
				XMin = min(v.X, XMin);
				XMax = max(v.X, XMax);
				YMin = min(v.Y, YMin);
				YMax = max(v.Y, YMax);
				ZMin = min(v.Z, ZMin);
				ZMax = max(v.Z, ZMax);
				vertices.push_back(v);
			}
			if (strcmp(token1, "vn") == 0) {
				CVec3Df N;

				sscanf(st2.NextTokenIE(" "), "%f", &N.X);
				sscanf(st2.NextTokenIE(" "), "%f", &N.Y);
				sscanf(st2.NextTokenIE(" "), "%f", &N.Z);
			}
		}
		if (strcmp(token1, "f") == 0) {
			CFace f;
			int *v = (int *)malloc(sizeof(int) * 1);
			int *tv = (int *)malloc(sizeof(int) * 1);
			int *n = (int *)malloc(sizeof(int) * 1);
			int cnt = 0;
			int sz = 1;
			token1 = st2.NextTokenIE(" ");
			while (token1 != NULL) {
				char *tokens2 = strdup(token1);
				char *token2;
				CStringTokenizer st3(tokens2);

				if (cnt == sz) {
					sz <<= 1;
					v = (int *)realloc(v, sizeof(int) * sz);
					tv = (int *)realloc(tv, sizeof(int) * sz);
					n = (int *)realloc(n, sizeof(int) * sz);
				}
				sscanf(st3.NextTokenIE("/"), "%d", &v[cnt]);
				token2 = st3.NextTokenAE("/");
				if (strlen(token2) > 0) sscanf(token2, "%d", &tv[cnt]);
				token2 = st3.NextTokenAE("/");
				if (token2 != NULL) sscanf(token2, "%d", &n[cnt]);
				token1 = st2.NextTokenIE(" ");
				++cnt;
			}
			for (int i = 0; i < cnt - 2; ++i) {
				f.P1 = v[0] - 1; f.P2 = v[i + 1] - 1; f.P3 = v[i + 2] - 1;
				faces.push_back(f);
			}
			free(v);
			free(tv);
			free(n);
		}
		line = st1.NextTokenIE("\n\r");
		free(tokens1);
	}
	free(buf);
	
	std::vector<CVec3Df> verticesV(vertices.begin(), vertices.end());
	std::vector<CFace> facesV(faces.begin(), faces.end());
	float *cdf = (float *)malloc(sizeof(float) * facesV.size());
	float A = 0.0f;
	for (int i = 0; i < facesV.size(); ++i) {
	    CFace f = facesV[i];
	    CVec3Df P1 = verticesV[f.P1];
	    CVec3Df P2 = verticesV[f.P2];
	    CVec3Df P3 = verticesV[f.P3];
	    CVec3Df U(P2.X - P1.X, P2.Y - P1.Y, P2.Z - P1.Z);
	    CVec3Df V(P3.X - P1.X, P3.Y - P1.Y, P3.Z - P1.Z);
	    CVec3Df N((U.Y * V.Z) - (V.Y * U.Z), (U.Z * V.X) - (V.Z * U.X), (U.X * V.Y) - (V.X * U.Y));
	    float A = 0.5f * sqrt((N.X * N.X) + (N.Y * N.Y) + (N.Z * N.Z));
	    if (i == 0) cdf[i] = A;
	    else
	        cdf[i] = cdf[i - 1] + A;
	}
	for (int i = 0; i < facesV.size(); ++i) cdf[i] /= cdf[facesV.size() - 1];
	cdf[facesV.size() - 1] = 1.0f;
	mat points(pointsNum, 3);
	for (int i = 0; i < pointsNum; ++i) {
	    int ind = InverseTransformSampling(cdf, facesV.size());
	    CVec3Df P = SampleTriUniform(verticesV, facesV[ind]);
	    points(i, 0) = P.X;
	    points(i, 1) = P.Y;
	    points(i, 2) = P.Z;
	}
	free(cdf);
	return points;
}