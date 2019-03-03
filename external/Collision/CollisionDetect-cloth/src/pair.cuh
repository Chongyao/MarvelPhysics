
typedef struct _g_pair {
	int2 *_dPairs;
	uint *_dIdx;
	int _offset;
	int _length;
	int _max_length;

	void init(int length) {
		uint dummy[] = { 0 };
		cutilSafeCall(cudaMalloc((void**)&_dIdx, 1 * sizeof(uint)));
		cutilSafeCall(cudaMemcpy(_dIdx, dummy, 1 * sizeof(uint), cudaMemcpyHostToDevice));

		_length = length;
		cutilSafeCall(cudaMalloc((void**)&_dPairs, length*sizeof(uint)));
		cutilSafeCall(cudaMemset(_dPairs, 0, length*sizeof(uint)));
		reportMemory("g_pair.init");

		_offset = 0;
		_max_length=0;
	}

	void append(_g_pair &pairs){
		uint len0[1];
		cutilSafeCall(cudaMemcpy(len0, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		uint len1[1];
		cutilSafeCall(cudaMemcpy(len1, pairs._dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		uint newlen[1];
		cutilSafeCall(cudaMemcpy(_dPairs + len0[0], pairs._dPairs, len1[0] * sizeof(int2), cudaMemcpyDeviceToDevice));
		newlen[0] = len0[0] + len1[0];
		cutilSafeCall(cudaMemcpy(_dIdx, newlen, 1 * sizeof(uint), cudaMemcpyHostToDevice));
	}

	void clear() {
		uint dummy[] = { 0 };
		cutilSafeCall(cudaMemcpy(_dIdx, dummy, 1 * sizeof(uint), cudaMemcpyHostToDevice));
		_offset = 0;
	}

	int getCollisions(bool self, struct _g_pair &rets,double *time);
	int getProximityConstraints(bool self, REAL mu, REAL mu_obs, REAL mrt, REAL mcs);
	int getImpacts(bool self, REAL mu, REAL mu_obs, _g_pair &vfPairs, _g_pair &eePairs, int &vfLen, int &eeLen);

	void destroy() {
		cudaFree(_dPairs);
		cudaFree(_dIdx);
	}

	uint length() {
		uint dummy[] = { 0 };
		cutilSafeCall(cudaMemcpy(dummy, _dIdx, 1 * sizeof(uint), cudaMemcpyDeviceToHost));
		if (dummy[0] > _max_length)
			_max_length = dummy[0];

		return dummy[0];
	}

	void setLength(uint len) {
		cutilSafeCall(cudaMemcpy(_dIdx, &len, 1 * sizeof(uint), cudaMemcpyHostToDevice));
	}

	int maxLength() const { return _max_length; }
} g_pair;
