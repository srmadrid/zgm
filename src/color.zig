const std = @import("std");
const zgm = @import("zgm.zig");

/// A color in the RGB color space.
pub const RGB = struct {
    v: @Vector(3, f32),

    /// Create a new color from the given red, green, and blue components. Each
    /// component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `rs`: The red component.
    /// - `gs`: The green component.
    /// - `bs`: The blue component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub inline fn init(rs: f32, gs: f32, bs: f32) RGB {
        return RGB{ .v = @Vector(3, f32){ rs, gs, bs } };
    }

    /// Create a new color from the given red, green, and blue components. Each
    /// component should be in the range [0, 255].
    ///
    /// **Parameters**:
    /// - `rs`: The red component.
    /// - `gs`: The green component.
    /// - `bs`: The blue component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub inline fn fromu8(rs: u8, gs: u8, bs: u8) RGB {
        return RGB{ .v = @Vector(3, f32){ @as(f32, @floatFromInt(rs)) / 255.0, @as(f32, @floatFromInt(gs)) / 255.0, @as(f32, @floatFromInt(bs)) / 255.0 } };
    }

    pub inline fn r(self: *const RGB) f32 {
        return self.v[0];
    }

    pub inline fn g(self: *const RGB) f32 {
        return self.v[1];
    }

    pub inline fn b(self: *const RGB) f32 {
        return self.v[2];
    }

    /// Convert the RGB color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    /// - `a`: The alpha component of the resulting RGBA color. Must be in the
    ///        range [0, 1].
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given RGB color.
    pub inline fn toRGBA(c: *const RGB, a: f32) RGBA {
        return RGBA{ .v = @Vector(4, f32){ c.v[0], c.v[1], c.v[2], a } };
    }

    /// Convert the RGB color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given RGB color.
    pub inline fn toHSL(c: *const RGB) HSL {
        const max = @max(@max(c.v[0], c.v[1]), c.v[2]);
        const min = @min(@min(c.v[0], c.v[1]), c.v[2]);
        const delta = max - min;

        var h: f32 = 0.0;
        var s: f32 = 0.0;
        const l: f32 = (max + min) / 2.0;

        if (delta != 0.0) {
            s = delta / (1.0 - @abs(2.0 * l - 1.0));

            if (max == c.v[0]) {
                h = (c.v[1] - c.v[2]) / delta;
            } else if (max == c.v[1]) {
                h = 2.0 + (c.v[2] - c.v[0]) / delta;
            } else {
                h = 4.0 + (c.v[0] - c.v[1]) / delta;
            }

            h *= 60.0;
            if (h < 0.0) {
                h += 360.0;
            }
        }

        return HSL{ .v = @Vector(3, f32){ h / 360.0, s, l } };
    }

    /// Convert the RGB color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    /// - `as`: The alpha component of the resulting HSLA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given RGB color.
    pub inline fn toHSLA(c: *const RGB, as: f32) HSLA {
        const hsl = c.toHSL();
        return HSLA{ .v = @Vector(4, f32){ hsl.v[0], hsl.v[1], hsl.v[2], as } };
    }

    /// Convert the RGB color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given RGB color.
    pub inline fn toHSV(c: *const RGB) HSV {
        const max = @max(@max(c.v[0], c.v[1]), c.v[2]);
        const min = @min(@min(c.v[0], c.v[1]), c.v[2]);
        const delta = max - min;

        var h: f32 = 0.0;
        var s: f32 = 0.0;
        const v: f32 = max;

        if (delta != 0.0) {
            s = delta / v;

            if (max == c.v[0]) {
                h = (c.v[1] - c.v[2]) / delta;
            } else if (max == c.v[1]) {
                h = 2.0 + (c.v[2] - c.v[0]) / delta;
            } else {
                h = 4.0 + (c.v[0] - c.v[1]) / delta;
            }

            h *= 60.0;
            if (h < 0.0) {
                h += 360.0;
            }
        }

        return HSV{ .vv = @Vector(3, f32){ h / 360.0, s, v } };
    }

    /// Convert the RGB color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    /// - `as`: The alpha component of the resulting HSVA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given RGB color.
    pub inline fn toHSVA(c: *const RGB, as: f32) HSVA {
        const hsv = c.toHSV();
        return HSVA{ .vv = @Vector(4, f32){ hsv.vv[0], hsv.vv[1], hsv.vv[2], as } };
    }

    /// Convert the RGB color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGB color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given RGB color.
    pub inline fn toCMYK(c: *const RGB) CMYK {
        const k = 1.0 - @max(@max(c.v[0], c.v[1]), c.v[2]);
        const cmy = @Vector(3, f32){ (1.0 - c.v[0] - k) / (1.0 - k), (1.0 - c.v[1] - k) / (1.0 - k), (1.0 - c.v[2] - k) / (1.0 - k) };
        return CMYK{ .v = @Vector(4, f32){ cmy[0], cmy[1], cmy[2], k } };
    }

    /// Compares two RGB colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const RGB, d: *const RGB) bool {
        return c.v[0] == d.v[0] and c.v[1] == d.v[1] and c.v[2] == d.v[2];
    }

    /// Compares two RGB colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const RGB, d: *const RGB, tolerance: f32) bool {
        return @abs(c.v[0] - d.v[0]) <= tolerance and @abs(c.v[1] - d.v[1]) <= tolerance and @abs(c.v[2] - d.v[2]) <= tolerance;
    }
};

/// A color in the RGBA color space.
pub const RGBA = struct {
    v: @Vector(4, f32),

    /// Create a new color from the given red, green, blue, and alpha components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `rs`: The red component.
    /// - `gs`: The green component.
    /// - `bs`: The blue component.
    /// - `as`: The alpha component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(rs: f32, gs: f32, bs: f32, as: f32) RGBA {
        return RGBA{ .v = @Vector(4, f32){ rs, gs, bs, as } };
    }

    /// Create a new color from the given red, green, blue, and alpha components.
    /// Each component should be in the range [0, 255].
    ///
    /// **Parameters**:
    /// - `rs`: The red component.
    /// - `gs`: The green component.
    /// - `bs`: The blue component.
    /// - `as`: The alpha component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn fromu8(rs: u8, gs: u8, bs: u8, as: u8) RGBA {
        return RGBA{ .v = @Vector(4, f32){ @as(f32, @floatFromInt(rs)) / 255.0, @as(f32, @floatFromInt(gs)) / 255.0, @as(f32, @floatFromInt(bs)) / 255.0, @as(f32, @floatFromInt(as)) / 255.0 } };
    }

    pub inline fn r(self: *const RGBA) f32 {
        return self.v[0];
    }

    pub inline fn g(self: *const RGBA) f32 {
        return self.v[1];
    }

    pub inline fn b(self: *const RGBA) f32 {
        return self.v[2];
    }

    pub inline fn a(self: *const RGBA) f32 {
        return self.v[3];
    }

    /// Convert the RGBA color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given RGBA color.
    pub inline fn toRGB(c: *const RGBA) RGB {
        return RGB{ .v = @Vector(3, f32){ c.v[0], c.v[1], c.v[2] } };
    }

    /// Convert the RGBA color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given RGBA color.
    pub inline fn toHSL(c: *const RGBA) HSL {
        return c.toRGB().toHSL();
    }

    /// Convert the RGBA color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given RGBA color.
    pub inline fn toHSLA(c: *const RGBA) HSLA {
        return c.toRGB().toHSLA(c.v[3]);
    }

    /// Convert the RGBA color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given RGBA color.
    pub inline fn toHSV(c: *const RGBA) HSV {
        return c.toRGB().toHSV();
    }

    /// Convert the RGBA color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given RGBA color.
    pub inline fn toHSVA(c: *const RGBA) HSVA {
        return c.toRGB().toHSVA(c.v[3]);
    }

    /// Convert the RGBA color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The RGBA color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given RGBA color.
    pub inline fn toCMYK(c: *const RGBA) CMYK {
        return c.toRGB().toCMYK();
    }

    /// Compares two RGBA colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const RGBA, d: *const RGBA) bool {
        return c.v[0] == d.v[0] and c.v[1] == d.v[1] and c.v[2] == d.v[2] and c.v[3] == d.v[3];
    }

    /// Compares two RGBA colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const RGBA, d: *const RGBA, tolerance: f32) bool {
        return @abs(c.v[0] - d.v[0]) <= tolerance and @abs(c.v[1] - d.v[1]) <= tolerance and @abs(c.v[2] - d.v[2]) <= tolerance and @abs(c.v[3] - d.v[3]) <= tolerance;
    }
};

/// A color in the HSL color space.
pub const HSL = struct {
    v: @Vector(3, f32),

    /// Create a new color from the given hue, saturation, and lightness components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `hs`: The hue component.
    /// - `ss`: The saturation component.
    /// - `ls`: The lightness component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(hs: f32, ss: f32, ls: f32) HSL {
        return HSL{ .v = @Vector(3, f32){ hs, ss, ls } };
    }

    pub inline fn h(self: *const HSL) f32 {
        return self.v[0];
    }

    pub inline fn s(self: *const HSL) f32 {
        return self.v[1];
    }

    pub inline fn l(self: *const HSL) f32 {
        return self.v[2];
    }

    /// Convert the HSL color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given HSL color.
    pub inline fn toRGB(c: *const HSL) RGB {
        if (c.v[1] == 0.0) {
            return RGB.init(c.v[2], c.v[2], c.v[2]);
        }

        const Hx6 = c.v[0] * 6.0;
        const C = (1.0 - @abs(2.0 * c.v[2] - 1.0)) * c.v[1];
        const X = C * (1.0 - @abs(@mod(Hx6, 2) - 1.0));
        const m = c.v[2] - C / 2.0;

        if (Hx6 < 1.0) {
            return RGB{ .v = @Vector(3, f32){ C + m, X + m, m } };
        } else if (Hx6 < 2.0) {
            return RGB{ .v = @Vector(3, f32){ X + m, C + m, m } };
        } else if (Hx6 < 3.0) {
            return RGB{ .v = @Vector(3, f32){ m, C + m, X + m } };
        } else if (Hx6 < 4.0) {
            return RGB{ .v = @Vector(3, f32){ m, X + m, C + m } };
        } else if (Hx6 < 5.0) {
            return RGB{ .v = @Vector(3, f32){ X + m, m, C + m } };
        } else {
            return RGB{ .v = @Vector(3, f32){ C + m, m, X + m } };
        }
    }

    /// Convert the HSL color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    /// - `as`: The alpha component of the resulting RGBA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given HSL color.
    pub inline fn toRGBA(c: *const HSL, as: f32) RGBA {
        return c.toRGB().toRGBA(as);
    }

    /// Convert the HSL color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    /// - `as`: The alpha component of the resulting HSLA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given HSL color.
    pub inline fn toHSLA(c: *const HSL, as: f32) HSLA {
        return HSLA{ .v = @Vector(4, f32){ c.v[0], c.v[1], c.v[2], as } };
    }

    /// Convert the HSL color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given HSL color.
    pub inline fn toHSV(c: *const HSL) HSV {
        const v = c.v[2] + c.v[1] * @min(c.v[2], 1.0 - c.v[2]);

        var ss: f32 = undefined;
        if (v == 0.0) {
            ss = 0.0;
        } else {
            ss = 2.0 * (v - c.v[2]) / v;
        }

        return HSV{ .vv = @Vector(3, f32){ c.v[0], ss, v } };
    }

    /// Convert the HSL color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    /// - `as`: The alpha component of the resulting HSVA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given HSL color.
    pub inline fn toHSVA(c: *const HSL, as: f32) HSVA {
        const hsv = c.toHSV();
        return HSVA{ .vv = @Vector(4, f32){ hsv.vv[0], hsv.vv[1], hsv.vv[2], as } };
    }

    /// Convert the HSL color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSL color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given HSL color.
    pub inline fn toCMYK(c: *const HSL) CMYK {
        return c.toRGB().toCMYK();
    }

    /// Compares two HSL colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const HSL, d: *const HSL) bool {
        return c.v[0] == d.v[0] and c.v[1] == d.v[1] and c.v[2] == d.v[2];
    }

    /// Compares two HSL colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const HSL, d: *const HSL, tolerance: f32) bool {
        return @abs(c.v[0] - d.v[0]) <= tolerance and @abs(c.v[1] - d.v[1]) <= tolerance and @abs(c.v[2] - d.v[2]) <= tolerance;
    }
};

/// A color in the HSLA color space.
pub const HSLA = struct {
    v: @Vector(4, f32),

    /// Create a new color from the given hue, saturation, lightness, and alpha components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `hs`: The hue component.
    /// - `ss`: The saturation component.
    /// - `ls`: The lightness component.
    /// - `as`: The alpha component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(hs: f32, ss: f32, ls: f32, as: f32) HSLA {
        return HSLA{ .v = @Vector(4, f32){ hs, ss, ls, as } };
    }

    pub inline fn h(self: *const HSLA) f32 {
        return self.v[0];
    }

    pub inline fn s(self: *const HSLA) f32 {
        return self.v[1];
    }

    pub inline fn l(self: *const HSLA) f32 {
        return self.v[2];
    }

    pub inline fn a(self: *const HSLA) f32 {
        return self.v[3];
    }

    /// Convert the HSLA color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given HSLA color.
    pub inline fn toRGB(c: *const HSLA) RGB {
        return c.toHSL().toRGB();
    }

    /// Convert the HSLA color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given HSLA color.
    pub inline fn toRGBA(c: *const HSLA) RGBA {
        return c.toHSL().toRGBA(c.v[3]);
    }

    /// Convert the HSLA color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given HSLA color.
    pub inline fn toHSL(c: *const HSLA) HSL {
        return HSL{ .v = @Vector(3, f32){ c.v[0], c.v[1], c.v[2] } };
    }

    /// Convert the HSLA color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given HSLA color.
    pub inline fn toHSV(c: *const HSLA) HSV {
        return c.toHSL().toHSV();
    }

    /// Convert the HSLA color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given HSLA color.
    pub inline fn toHSVA(c: *const HSLA) HSVA {
        return c.toHSL().toHSVA(c.v[3]);
    }

    /// Convert the HSLA color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSLA color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given HSLA color.
    pub inline fn toCMYK(c: *const HSLA) CMYK {
        return c.toRGB().toCMYK();
    }

    /// Compares two HSLA colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const HSLA, d: *const HSLA) bool {
        return c.v[0] == d.v[0] and c.v[1] == d.v[1] and c.v[2] == d.v[2] and c.v[3] == d.v[3];
    }

    /// Compares two HSLA colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const HSLA, d: *const HSLA, tolerance: f32) bool {
        return @abs(c.v[0] - d.v[0]) <= tolerance and @abs(c.v[1] - d.v[1]) <= tolerance and @abs(c.v[2] - d.v[2]) <= tolerance and @abs(c.v[3] - d.v[3]) <= tolerance;
    }
};

/// A color in the HSV color space.
pub const HSV = struct {
    vv: @Vector(3, f32),

    /// Create a new color from the given hue, saturation, and value components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `hs`: The hue component.
    /// - `ss`: The saturation component.
    /// - `vs`: The value component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(hs: f32, ss: f32, vs: f32) HSV {
        return HSV{ .vv = @Vector(3, f32){ hs, ss, vs } };
    }

    pub inline fn h(self: *const HSV) f32 {
        return self.vv[0];
    }

    pub inline fn s(self: *const HSV) f32 {
        return self.vv[1];
    }

    pub inline fn v(self: *const HSV) f32 {
        return self.vv[2];
    }

    /// Convert the HSV color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given HSV color.
    pub inline fn toRGB(c: *const HSV) RGB {
        if (c.vv[1] == 0.0) {
            return RGB.init(c.vv[2], c.vv[2], c.vv[2]);
        }

        const Hx6 = c.vv[0] * 6.0;
        const C = c.vv[2] * c.vv[1];
        const X = C * (1.0 - @abs(@mod(Hx6, 2) - 1.0));
        const m = c.vv[2] - C;

        if (Hx6 < 1.0) {
            return RGB{ .v = @Vector(3, f32){ C + m, X + m, m } };
        } else if (Hx6 < 2.0) {
            return RGB{ .v = @Vector(3, f32){ X + m, C + m, m } };
        } else if (Hx6 < 3.0) {
            return RGB{ .v = @Vector(3, f32){ m, C + m, X + m } };
        } else if (Hx6 < 4.0) {
            return RGB{ .v = @Vector(3, f32){ m, X + m, C + m } };
        } else if (Hx6 < 5.0) {
            return RGB{ .v = @Vector(3, f32){ X + m, m, C + m } };
        } else {
            return RGB{ .v = @Vector(3, f32){ C + m, m, X + m } };
        }
    }

    /// Convert the HSV color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    /// - `as`: The alpha component of the resulting RGBA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given HSV color.
    pub inline fn toRGBA(c: *const HSV, as: f32) RGBA {
        return c.toRGB().toRGBA(as);
    }

    /// Convert the HSV color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given HSV color.
    pub inline fn toHSL(c: *const HSV) HSL {
        const l = c.vv[2] * (1.0 - c.vv[1] / 2.0);

        var ss: f32 = undefined;
        if (l == 0.0 or l == 1.0) {
            ss = 0.0;
        } else {
            ss = (c.vv[2] - l) / @min(l, 1.0 - l);
        }

        return HSL{ .v = @Vector(3, f32){ c.vv[0], ss, l } };
    }

    /// Convert the HSV color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    /// - `as`: The alpha component of the resulting HSLA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given HSV color.
    pub inline fn toHSLA(c: *const HSV, as: f32) HSLA {
        return c.toHSL().toHSLA(as);
    }

    /// Convert the HSV color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    /// - `as`: The alpha component of the resulting HSVA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given HSV color.
    pub inline fn toHSVA(c: *const HSV, as: f32) HSVA {
        return HSVA{ .vv = @Vector(4, f32){ c.vv[0], c.vv[1], c.vv[2], as } };
    }

    /// Convert the HSV color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSV color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given HSV color.
    pub inline fn toCMYK(c: *const HSV) CMYK {
        return c.toRGB().toCMYK();
    }

    /// Compares two HSV colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const HSV, d: *const HSV) bool {
        return c.vv[0] == d.vv[0] and c.vv[1] == d.vv[1] and c.vv[2] == d.vv[2];
    }

    /// Compares two HSV colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const HSV, d: *const HSV, tolerance: f32) bool {
        return @abs(c.vv[0] - d.vv[0]) <= tolerance and @abs(c.vv[1] - d.vv[1]) <= tolerance and @abs(c.vv[2] - d.vv[2]) <= tolerance;
    }
};

/// A color in the HSVA color space.
pub const HSVA = struct {
    vv: @Vector(4, f32),

    /// Create a new color from the given hue, saturation, value, and alpha components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `hs`: The hue component.
    /// - `ss`: The saturation component.
    /// - `vs`: The value component.
    /// - `as`: The alpha component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(hs: f32, ss: f32, vs: f32, as: f32) HSVA {
        return HSVA{ .vv = @Vector(4, f32){ hs, ss, vs, as } };
    }

    pub inline fn h(self: *const HSVA) f32 {
        return self.vv[0];
    }

    pub inline fn s(self: *const HSVA) f32 {
        return self.vv[1];
    }

    pub inline fn v(self: *const HSVA) f32 {
        return self.vv[2];
    }

    pub inline fn a(self: *const HSVA) f32 {
        return self.vv[3];
    }

    /// Convert the HSVA color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given HSVA color.
    pub inline fn toRGB(c: *const HSVA) RGB {
        return c.toHSV().toRGB();
    }

    /// Convert the HSVA color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given HSVA color.
    pub inline fn toRGBA(c: *const HSVA) RGBA {
        return c.toHSV().toRGBA(c.vv[3]);
    }

    /// Convert the HSVA color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given HSVA color.
    pub inline fn toHSL(c: *const HSVA) HSL {
        return c.toHSV().toHSL();
    }

    /// Convert the HSVA color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given HSVA color.
    pub inline fn toHSLA(c: *const HSVA) HSLA {
        return c.toHSV().toHSLA(c.vv[3]);
    }

    /// Convert the HSVA color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given HSVA color.
    pub inline fn toHSV(c: *const HSVA) HSV {
        return HSV{ .vv = @Vector(3, f32){ c.vv[0], c.vv[1], c.vv[2] } };
    }

    /// Convert the HSVA color to the CMYK color space.
    ///
    /// **Parameters**:
    /// - `c`: The HSVA color to convert.
    ///
    /// **Returns**:
    /// - The CMYK color equivalent to the given HSVA color.
    pub inline fn toCMYK(c: *const HSVA) CMYK {
        return c.toRGB().toCMYK();
    }

    /// Compares two HSVA colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(c: *const HSVA, d: *const HSVA) bool {
        return c.vv[0] == d.vv[0] and c.vv[1] == d.vv[1] and c.vv[2] == d.vv[2] and c.vv[3] == d.vv[3];
    }

    /// Compares two HSVA colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `c`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(c: *const HSVA, d: *const HSVA, tolerance: f32) bool {
        return @abs(c.vv[0] - d.vv[0]) <= tolerance and @abs(c.vv[1] - d.vv[1]) <= tolerance and @abs(c.vv[2] - d.vv[2]) <= tolerance and @abs(c.vv[3] - d.vv[3]) <= tolerance;
    }
};

/// A color in the CMYK color space.
pub const CMYK = struct {
    v: @Vector(4, f32),

    /// Create a new color from the given cyan, magenta, yellow, and key components.
    /// Each component should be in the range [0, 1].
    ///
    /// **Parameters**:
    /// - `cs`: The cyan component.
    /// - `ms`: The magenta component.
    /// - `ys`: The yellow component.
    /// - `ks`: The key component.
    ///
    /// **Returns**:
    /// - A new color with the given components.
    pub fn init(cs: f32, ms: f32, ys: f32, ks: f32) CMYK {
        return CMYK{ .v = @Vector(4, f32){ cs, ms, ys, ks } };
    }

    pub inline fn c(self: *const CMYK) f32 {
        return self.v[0];
    }

    pub inline fn m(self: *const CMYK) f32 {
        return self.v[1];
    }

    pub inline fn y(self: *const CMYK) f32 {
        return self.v[2];
    }

    pub inline fn k(self: *const CMYK) f32 {
        return self.v[3];
    }

    /// Convert the CMYK color to the RGB color space.
    ///
    /// **Parameters**:
    /// - `c`: The CMYK color to convert.
    ///
    /// **Returns**:
    /// - The RGB color equivalent to the given CMYK color.
    pub inline fn toRGB(cc: *const CMYK) RGB {
        return RGB{ .v = @Vector(3, f32){ (1.0 - cc.v[0]) * (1.0 - cc.v[3]), (1.0 - cc.v[1]) * (1.0 - cc.v[3]), (1.0 - cc.v[2]) * (1.0 - cc.v[3]) } };
    }

    /// Convert the CMYK color to the RGBA color space.
    ///
    /// **Parameters**:
    /// - `cc`: The CMYK color to convert.
    /// - `as`: The alpha component of the resulting RGBA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The RGBA color equivalent to the given CMYK color.
    pub inline fn toRGBA(cc: *const CMYK, as: f32) RGBA {
        return cc.toRGB().toRGBA(as);
    }

    /// Convert the CMYK color to the HSL color space.
    ///
    /// **Parameters**:
    /// - `cc`: The CMYK color to convert.
    ///
    /// **Returns**:
    /// - The HSL color equivalent to the given CMYK color.
    pub inline fn toHSL(cc: *const CMYK) HSL {
        return cc.toRGB().toHSL();
    }

    /// Convert the CMYK color to the HSLA color space.
    ///
    /// **Parameters**:
    /// - `cc`: The CMYK color to convert.
    /// - `as`: The alpha component of the resulting HSLA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSLA color equivalent to the given CMYK color.
    pub inline fn toHSLA(cc: *const CMYK, as: f32) HSLA {
        return cc.toRGB().toHSLA(as);
    }

    /// Convert the CMYK color to the HSV color space.
    ///
    /// **Parameters**:
    /// - `cc`: The CMYK color to convert.
    ///
    /// **Returns**:
    /// - The HSV color equivalent to the given CMYK color.
    pub inline fn toHSV(cc: *const CMYK) HSV {
        return cc.toRGB().toHSV();
    }

    /// Convert the CMYK color to the HSVA color space.
    ///
    /// **Parameters**:
    /// - `cc`: The CMYK color to convert.
    /// - `as`: The alpha component of the resulting HSVA color. Must be in the
    ///         range [0, 1].
    ///
    /// **Returns**:
    /// - The HSVA color equivalent to the given CMYK color.
    pub inline fn toHSVA(cc: *const CMYK, as: f32) HSVA {
        return cc.toRGB().toHSVA(as);
    }

    /// Compares two CMYK colors and returns true if they are equal.
    ///
    /// **Parameters**:
    /// - `cc`: The first color.
    /// - `d`: The second color.
    ///
    /// **Returns**:
    /// - `true` if the colors are equal, `false` otherwise.
    pub inline fn equal(cc: *const CMYK, d: *const CMYK) bool {
        return cc.v[0] == d.v[0] and cc.v[1] == d.v[1] and cc.v[2] == d.v[2] and cc.v[3] == d.v[3];
    }

    /// Compares two CMYK colors and returns true if they are approximately
    /// equal.
    ///
    /// **Parameters**:
    /// - `cc`: The first color.
    /// - `d`: The second color.
    /// - `tolerance`: The maximum difference between the components.
    ///
    /// **Returns**:
    /// - `true` if the colors are approximately equal, `false` otherwise.
    pub inline fn approxEqual(cc: *const CMYK, d: *const CMYK, tolerance: f32) bool {
        return @abs(cc.v[0] - d.v[0]) <= tolerance and @abs(cc.v[1] - d.v[1]) <= tolerance and @abs(cc.v[2] - d.v[2]) <= tolerance and @abs(cc.v[3] - d.v[3]) <= tolerance;
    }
};

test "RGB.toRGBA" {
    const rgb = RGB.init(0.5, 0.25, 0.75);
    const r = rgb.toRGBA(0.5);

    try std.testing.expect(r.approxEqual(&RGBA.init(0.5, 0.25, 0.75, 0.5), 0.0001));
}

test "RGB.toHSL" {
    const rgb = RGB.fromu8(255, 128, 64);
    const r = rgb.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.0001));
}

test "RGB.toHSLA" {
    const rgb = RGB.fromu8(255, 128, 64);
    const r = rgb.toHSLA(0.5);

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5), 0.0001));
}

test "RGB.toHSV" {
    const rgb = RGB.fromu8(255, 128, 64);
    const r = rgb.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.0001));
}

test "RGB.toHSVA" {
    const rgb = RGB.fromu8(255, 128, 64);
    const r = rgb.toHSVA(0.5);

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5), 0.0001));
}

test "RGB.toCMYK" {
    const rgb = RGB.fromu8(255, 128, 64);
    const r = rgb.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.498, 0.749, 0.0), 0.0001));
}

test "RGBA.toRGB" {
    const rgba = RGBA.init(0.5, 0.25, 0.75, 0.5);
    const r = rgba.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.init(0.5, 0.25, 0.75), 0.0001));
}

test "RGBA.toHSL" {
    const rgba = RGBA.fromu8(255, 128, 64, 128);
    const r = rgba.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.0001));
}

test "RGBA.toHSLA" {
    const rgba = RGBA.fromu8(255, 128, 64, 128);
    const r = rgba.toHSLA();

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5020), 0.0001));
}

test "RGBA.toHSV" {
    const rgba = RGBA.fromu8(255, 128, 64, 128);
    const r = rgba.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.0001));
}

test "RGBA.toHSVA" {
    const rgba = RGBA.fromu8(255, 128, 64, 128);
    const r = rgba.toHSVA();

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5020), 0.0001));
}

test "RGBA.toCMYK" {
    const rgba = RGBA.fromu8(255, 128, 64, 128);
    const r = rgba.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.498, 0.749, 0.0), 0.0001));
}

test "HSL.toRGB" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.fromu8(255, 128, 64), 0.001));
}

test "HSL.toRGBA" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toRGBA(0.5020);

    try std.testing.expect(r.approxEqual(&RGBA.fromu8(255, 128, 64, 128), 0.001));
}

test "HSL.toHSLA" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toHSLA(0.5020);

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5020), 0.001));
}

test "HSL.toHSV" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.001));
}

test "HSL.toHSVA" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toHSVA(0.5020);

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5020), 0.001));
}

test "HSL.toCMYK" {
    const hsl = HSL.init(0.0558, 1.0, 0.6255);
    const r = hsl.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.4981, 0.749, 0.0), 0.001));
}

test "HSLA.toRGB" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.fromu8(255, 128, 64), 0.001));
}

test "HSLA.toRGBA" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toRGBA();

    try std.testing.expect(r.approxEqual(&RGBA.fromu8(255, 128, 64, 128), 0.001));
}

test "HSLA.toHSL" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.001));
}

test "HSLA.toHSV" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.001));
}

test "HSLA.toHSVA" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toHSVA();

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5020), 0.001));
}

test "HSLA.toCMYK" {
    const hsla = HSLA.init(0.0558, 1.0, 0.6255, 0.5020);
    const r = hsla.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.4981, 0.749, 0.0), 0.001));
}

test "HSV.toRGB" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.fromu8(255, 128, 64), 0.001));
}

test "HSV.toRGBA" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toRGBA(0.5020);

    try std.testing.expect(r.approxEqual(&RGBA.fromu8(255, 128, 64, 128), 0.001));
}

test "HSV.toHSL" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.001));
}

test "HSV.toHSLA" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toHSLA(0.5020);

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5020), 0.001));
}

test "HSV.toHSVA" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toHSVA(0.5020);

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5020), 0.001));
}

test "HSV.toCMYK" {
    const hsv = HSV.init(0.0558, 0.749, 1.0);
    const r = hsv.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.4981, 0.749, 0.0), 0.001));
}

test "HSVA.toRGB" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.fromu8(255, 128, 64), 0.001));
}

test "HSVA.toRGBA" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toRGBA();

    try std.testing.expect(r.approxEqual(&RGBA.fromu8(255, 128, 64, 128), 0.001));
}

test "HSVA.toHSL" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.001));
}

test "HSVA.toHSLA" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toHSLA();

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5020), 0.001));
}

test "HSVA.toHSV" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.001));
}

test "HSVA.toCMYK" {
    const hsva = HSVA.init(0.0558, 0.749, 1.0, 0.5020);
    const r = hsva.toCMYK();

    try std.testing.expect(r.approxEqual(&CMYK.init(0.0, 0.4981, 0.749, 0.0), 0.001));
}

test "CMYK.toRGB" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toRGB();

    try std.testing.expect(r.approxEqual(&RGB.fromu8(255, 128, 64), 0.001));
}

test "CMYK.toRGBA" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toRGBA(0.5020);

    try std.testing.expect(r.approxEqual(&RGBA.fromu8(255, 128, 64, 128), 0.001));
}

test "CMYK.toHSL" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toHSL();

    try std.testing.expect(r.approxEqual(&HSL.init(0.0558, 1.0, 0.6255), 0.001));
}

test "CMYK.toHSLA" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toHSLA(0.5020);

    try std.testing.expect(r.approxEqual(&HSLA.init(0.0558, 1.0, 0.6255, 0.5020), 0.001));
}

test "CMYK.toHSV" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toHSV();

    try std.testing.expect(r.approxEqual(&HSV.init(0.0558, 0.749, 1.0), 0.001));
}

test "CMYK.toHSVA" {
    const cmyk = CMYK.init(0.0, 0.4981, 0.749, 0.0);
    const r = cmyk.toHSVA(0.5020);

    try std.testing.expect(r.approxEqual(&HSVA.init(0.0558, 0.749, 1.0, 0.5020), 0.001));
}
