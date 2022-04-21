use criterion::BenchmarkId;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kmers::naive_impl;
use random_string::generate;

pub fn compute_rc_naive(b: &[u8]) -> u64 {
  let K = 31;
  let mut sum = 0u64;
  for i in 0..(b.len()-K) {
	let k = naive_impl::Kmer::from(&b[i..i+K]);
	//k.to_reverse_complement();
	sum += k.into_u64();
  }
  return sum;
}

pub fn compute_rc(b: &[u8]) -> u64 {
	const K:usize = 31;

	let mut sum = 0u64;
	let enc = kmers::encoding::xor10::Xor10;
	for i in 0..(b.len()-K) {
	      let k = kmers::kmer::Kmer::<u64, K, { kmers::kmer::word_for_k::<u64, K>() }>::new(&b[i..i+K], &enc);
	      sum += k.num_bytes() as u64;
	}
	return sum;
}

pub fn criterion_benchmark(c: &mut Criterion) {

    let charset = "ACGT";
    let input = generate(1_000_000, charset);
    let bytes = input.as_bytes();
    c.bench_with_input(BenchmarkId::new("construct naive independent", "1M"), &bytes, |b, &s| {
        b.iter(|| compute_rc_naive(&s));
    });

    c.bench_with_input(BenchmarkId::new("construct kme.rs independent", "1M"), &bytes, |b, &s| {
        b.iter(|| compute_rc(&s));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);