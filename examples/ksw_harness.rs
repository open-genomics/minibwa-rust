//! ksw2 对照工具：与 C 版 harness.c 输出格式一致。
use minibwa_rust::ksw2::{gen_nt4_mat, ksw_extd2, Extz};

fn encode(s: &str) -> Vec<u8> {
    s.bytes()
        .map(|b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let q: i8 = args[1].parse().unwrap();
    let e: i8 = args[2].parse().unwrap();
    let q2: i8 = args[3].parse().unwrap();
    let e2: i8 = args[4].parse().unwrap();
    let w: i32 = args[5].parse().unwrap();
    let zdrop: i32 = args[6].parse().unwrap();
    let end_bonus: i32 = args[7].parse().unwrap();
    let fstr = &args[8];
    let flag: u32 = if fstr.starts_with("0x") {
        u32::from_str_radix(&fstr[2..], 16).unwrap()
    } else {
        fstr.parse().unwrap()
    };
    let a: i8 = args[9].parse().unwrap();
    let b: i8 = args[10].parse().unwrap();
    let query = encode(&args[11]);
    let target = encode(&args[12]);
    let mut mat = [0i8; 25];
    gen_nt4_mat(&mut mat, a, if b < 0 { -b } else { b }, 0, 1, 0);
    let mut ez = Extz::default();
    ksw_extd2(
        query.len() as i32, &query, target.len() as i32, &target, 5, &mat,
        q, e, q2, e2, w, zdrop, end_bonus, flag, &mut ez,
    );
    print!(
        "score={} max={} max_q={} max_t={} mqe={} mqe_t={} zdropped={} reach_end={} cigar=",
        ez.score, ez.max, ez.max_q, ez.max_t, ez.mqe, ez.mqe_t,
        i32::from(ez.zdropped), i32::from(ez.reach_end)
    );
    for &c in &ez.cigar {
        print!("{}{}", c >> 4, "MIDN".as_bytes()[(c & 0xf) as usize] as char);
    }
    println!();
}
