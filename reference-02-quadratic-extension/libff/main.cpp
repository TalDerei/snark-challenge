#include <cstdio>
#include <vector>

#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_init.hpp>

using namespace libff;
using namespace std;

Fq<mnt4753_pp> read_mnt4_fq(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

void write_mnt4_fq(FILE* output, Fq<mnt4753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

Fqe<mnt4753_pp> read_mnt4_fq2(FILE* input) {
  Fq<mnt4753_pp> c0 = read_mnt4_fq(input);
  Fq<mnt4753_pp> c1 = read_mnt4_fq(input);
  return Fqe<mnt4753_pp>(c0, c1);
}

void write_mnt4_fq2(FILE* output, Fqe<mnt4753_pp> x) {
  write_mnt4_fq(output, x.c0);
  write_mnt4_fq(output, x.c1);
}

Fq<mnt4753_pp> read_mnt4_fq_numeral(FILE* input) {
  Fq<mnt4753_pp> x;
  // cout << "mnt4753_q_limbs is: " << mnt4753_q_limbs << endl;
  // cout << "sizeof(mp_size_t) is: " << sizeof(mp_size_t) << endl;
  // cout << "libff::mnt4753_q_limbs * sizeof(mp_size_t) is: " << libff::mnt4753_q_limbs * sizeof(mp_size_t) << endl;

  // 12 64-bit limbs = 768-bits (96 bytes) into x
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  // for (int i = 0; i < 12; i++) {
  //   cout << x.mont_repr.data[i] << endl;
  // }
  // exit(0);
  auto b = Fq<mnt4753_pp>(x.mont_repr);
  return b;
}

void write_mnt4_fq_numeral(FILE* output, Fq<mnt4753_pp> x) {
  auto out_numeral = x.as_bigint();
  fwrite((void *) out_numeral.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

Fqe<mnt4753_pp> read_mnt4_fq2_numeral(FILE* input) {
  Fq<mnt4753_pp> c0 = read_mnt4_fq_numeral(input);
  Fq<mnt4753_pp> c1 = read_mnt4_fq_numeral(input);
  // Create a pair c0,c1 which is how field elements are represented in Fqe finite field
  // Each element in the field (i.e. the pair (c0, c1) is of size 192 bytes)
  return Fqe<mnt4753_pp>(c0, c1);
}

void write_mnt4_fq2_numeral(FILE* output, Fqe<mnt4753_pp> x) {
  write_mnt4_fq_numeral(output, x.c0);
  write_mnt4_fq_numeral(output, x.c1);
}

int main(int argc, char *argv[])
{
    // argv should be
    // { "main", "compute" or "compute-numeral", inputs, outputs }

    mnt4753_pp::init_public_params();

    auto is_numeral = strcmp(argv[1], "compute-numeral") == 0;
    auto inputs = fopen(argv[2], "r");
    auto outputs = fopen(argv[3], "w");

    auto read_mnt4 = read_mnt4_fq;
    auto write_mnt4 = write_mnt4_fq;
    auto read_mnt4_q2 = read_mnt4_fq2;
    auto write_mnt4_q2 = write_mnt4_fq2;
    if (is_numeral) {
      read_mnt4 = read_mnt4_fq_numeral;
      write_mnt4 = write_mnt4_fq_numeral;
      read_mnt4_q2 = read_mnt4_fq2_numeral;
      write_mnt4_q2 = write_mnt4_fq2_numeral;
    }

    while (true) {
      size_t n;
      size_t elts_read = fread((void *) &n, sizeof(size_t), 1, inputs);

      if (elts_read == 0) { break; }

      // Define vectors x and y of q finite field elements from (Fq^2) mnt4753_pp, where each element 
      // is a pair (a0, a1) from Fq.
      // x be will 1024 pairs of 196 bytes each, where each pair can be further deconstructed to 12 8 byte limbs
      std::vector<Fqe<mnt4753_pp>> x;
      for (size_t i = 0; i < n; ++i) {
        x.emplace_back(read_mnt4_q2(inputs));
      }

      std::vector<Fqe<mnt4753_pp>> y;
      for (size_t i = 0; i < n; ++i) {
        y.emplace_back(read_mnt4_q2(inputs));
      }

      // The actual code for doing Fq2 multiplication lives in libff/algebra/fields/fp2.tcc
      for (size_t i = 0; i < n; ++i) {
        write_mnt4_q2(outputs, x[i] * y[i]);
      }
    }
    fclose(outputs);

    return 0;
}
