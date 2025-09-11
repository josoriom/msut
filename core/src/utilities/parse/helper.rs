#[inline]
pub fn set_u32_at(buf: &mut [u8], pos: usize, val: u32) {
    buf[pos..pos + 4].copy_from_slice(&val.to_le_bytes());
}
#[inline]
pub fn set_u64_at(buf: &mut [u8], pos: usize, val: u64) {
    buf[pos..pos + 8].copy_from_slice(&val.to_le_bytes());
}
#[inline]
pub fn set_f64_at(buf: &mut [u8], pos: usize, val: f64) {
    buf[pos..pos + 8].copy_from_slice(&val.to_le_bytes());
}
#[inline]
pub fn ensure_cap(out: &mut Vec<u8>, need: usize) {
    if need > out.len() {
        let mut cap = out.len().max(64);
        while cap < need {
            cap *= 2;
        }
        out.resize(cap, 0);
    }
}
#[inline]
pub unsafe fn write_f64_le(out: &mut Vec<u8>, pos: &mut usize, vals: &[f64]) -> (u64, u32) {
    *pos = (*pos + 7) & !7;
    let off = *pos as u64;
    let nbytes = vals.len() * 8;
    ensure_cap(out, *pos + nbytes);
    unsafe {
        std::ptr::copy_nonoverlapping(
            vals.as_ptr() as *const u8,
            out.as_mut_ptr().add(*pos),
            nbytes,
        )
    };
    *pos += nbytes;
    (off, vals.len() as u32)
}
#[inline]
pub unsafe fn write_f32_le(out: &mut Vec<u8>, pos: &mut usize, vals: &[f32]) -> (u64, u32) {
    *pos = (*pos + 7) & !7;
    let off = *pos as u64;
    let nbytes = vals.len() * 4;
    ensure_cap(out, *pos + nbytes);
    unsafe {
        std::ptr::copy_nonoverlapping(
            vals.as_ptr() as *const u8,
            out.as_mut_ptr().add(*pos),
            nbytes,
        )
    };
    *pos += nbytes;
    (off, vals.len() as u32)
}
#[inline]
pub unsafe fn write_f64_at(buf: &mut [u8], off: usize, vals: &[f64]) {
    unsafe {
        std::ptr::copy_nonoverlapping(
            vals.as_ptr() as *const u8,
            buf.as_mut_ptr().add(off),
            vals.len() * 8,
        )
    };
}
#[inline]
pub unsafe fn write_f32_at(buf: &mut [u8], off: usize, vals: &[f32]) {
    unsafe {
        std::ptr::copy_nonoverlapping(
            vals.as_ptr() as *const u8,
            buf.as_mut_ptr().add(off),
            vals.len() * 4,
        )
    };
}
#[inline]
pub fn rd_u32(b: &[u8], p: usize) -> Result<u32, String> {
    if p + 4 > b.len() {
        return Err("u32 OOB".into());
    }
    Ok(u32::from_le_bytes(b[p..p + 4].try_into().unwrap()))
}
#[inline]
pub fn rd_u64(b: &[u8], p: usize) -> Result<u64, String> {
    if p + 8 > b.len() {
        return Err("u64 OOB".into());
    }
    Ok(u64::from_le_bytes(b[p..p + 8].try_into().unwrap()))
}
#[inline]
pub fn rd_f64(b: &[u8], p: usize) -> Result<f64, String> {
    if p + 8 > b.len() {
        return Err("f64 OOB".into());
    }
    Ok(f64::from_le_bytes(b[p..p + 8].try_into().unwrap()))
}
#[inline]
pub fn read_array_as_f64(
    buf: &[u8],
    off: u64,
    len: u32,
    fmt: u8,
) -> Result<Option<Vec<f64>>, String> {
    if off == 0 || len == 0 {
        return Ok(None);
    }
    let off = off as usize;
    let n = len as usize;
    match fmt {
        2 => {
            let need = n.checked_mul(8).ok_or("len overflow")?;
            if off + need > buf.len() {
                return Err("f64 array OOB".into());
            }
            let mut out = Vec::with_capacity(n);
            let mut p = off;
            for _ in 0..n {
                out.push(f64::from_le_bytes(buf[p..p + 8].try_into().unwrap()));
                p += 8;
            }
            Ok(Some(out))
        }
        1 => {
            let need = n.checked_mul(4).ok_or("len overflow")?;
            if off + need > buf.len() {
                return Err("f32 array OOB".into());
            }
            let mut out = Vec::with_capacity(n);
            let mut p = off;
            for _ in 0..n {
                let v = f32::from_le_bytes(buf[p..p + 4].try_into().unwrap());
                out.push(v as f64);
                p += 4;
            }
            Ok(Some(out))
        }
        _ => Err("unknown fmt".into()),
    }
}
#[inline]
pub fn read_array_as_f32(
    buf: &[u8],
    off: u64,
    len: u32,
    fmt: u8,
) -> Result<Option<Vec<f32>>, String> {
    if off == 0 || len == 0 {
        return Ok(None);
    }
    let off = off as usize;
    let n = len as usize;
    #[inline]
    fn f64_to_f32_lossy(v: f64) -> f32 {
        if !v.is_finite() {
            return v as f32;
        }
        if v > f32::MAX as f64 {
            f32::MAX
        } else if v < -(f32::MAX as f64) {
            -f32::MAX
        } else {
            v as f32
        }
    }
    match fmt {
        1 => {
            let need = n.checked_mul(4).ok_or("len overflow")?;
            if off + need > buf.len() {
                return Err("f32 array OOB".into());
            }
            #[cfg(target_endian = "little")]
            unsafe {
                let mut out: Vec<f32> = Vec::with_capacity(n);
                out.set_len(n);
                std::ptr::copy_nonoverlapping(
                    buf.as_ptr().add(off),
                    out.as_mut_ptr() as *mut u8,
                    need,
                );
                return Ok(Some(out));
            }
            #[cfg(not(target_endian = "little"))]
            {
                let mut out = Vec::with_capacity(n);
                let mut p = off;
                for _ in 0..n {
                    out.push(f32::from_le_bytes(buf[p..p + 4].try_into().unwrap()));
                    p += 4;
                }
                Ok(Some(out))
            }
        }
        2 => {
            let need = n.checked_mul(8).ok_or("len overflow")?;
            if off + need > buf.len() {
                return Err("f64 array OOB".into());
            }
            let mut out = Vec::with_capacity(n);
            let mut p = off;
            for _ in 0..n {
                let v = f64::from_le_bytes(buf[p..p + 8].try_into().unwrap());
                out.push(f64_to_f32_lossy(v));
                p += 8;
            }
            Ok(Some(out))
        }
        _ => Err("unknown fmt".into()),
    }
}
