// Safety-layer curve construction. All candidates are verified by the caller.
namespace {

struct CornerPoint {
    V3 pos;
    V3 gun;
    double speed;
    int spray;
    int kind;
    int source;
};

bool sameProcess(const CornerPoint& a, const CornerPoint& b) {
    return a.spray == b.spray && a.kind == b.kind && a.speed == b.speed;
}

V3 blendAxis(V3 a, V3 b, double t) {
    const double cosine = std::clamp(dot(a, b), -1.0, 1.0);
    if (cosine > 0.999999) return unit(a * (1 - t) + b * t);
    V3 side = b - a * cosine;
    if (norm(side) < 1e-9)
        side = cross(a, std::fabs(a.z) < 0.9 ? V3(0, 0, 1) : V3(1, 0, 0));
    const double angle = t * std::acos(cosine);
    return unit(a * std::cos(angle) + unit(side) * std::sin(angle));
}

bool sampleAlongCurve(const std::vector<CornerPoint>& dense, double step,
                      std::vector<CornerPoint>& output) {
    output.clear();
    output.push_back(dense.front());
    size_t pending = 1;
    auto followsCurve = [&](V3 end, size_t limit) {
        const V3 start = output.back().pos;
        const V3 chord = end - start;
        const double squared = dot(chord, chord);
        for (; pending < limit; ++pending) {
            const V3 v = dense[pending].pos - start;
            const double fraction = squared > 1e-12 ? std::clamp(dot(v, chord) / squared, 0.0, 1.0) : 0.0;
            if (norm(v - chord * fraction) > 0.75) return false;
        }
        return true;
    };
    double distance = 0, next = step;
    for (size_t i = 1; i < dense.size(); ++i) {
        const auto& a = dense[i - 1];
        const auto& b = dense[i];
        const V3 start = a.pos, delta = b.pos - start;
        const double length = norm(delta);
        if (length <= 1e-10) continue;
        while (next <= distance + length - (i + 1 == dense.size() ? 1e-8 : 0.0)) {
            const double f = std::clamp((next - distance) / length, 0.0, 1.0);
            const V3 p = start + delta * f;
            if (!followsCurve(p, i)) return false;
            const V3 gun = blendAxis(a.gun, b.gun, f);
            output.push_back({p, gun, a.speed, a.spray, a.kind, f < 0.5 ? a.source : b.source});
            next += step;
        }
        distance += length;
    }
    if (!followsCurve(dense.back().pos, dense.size() - 1)) return false;
    output.push_back(dense.back());
    return true;
}

// Triangular kernel integrated against the input polyline in unsigned arc length.
// Its first two derivatives give continuous tangent and curvature across a corner.
class IntegratedPath {
public:
    IntegratedPath(const std::vector<CornerPoint>& path, const std::vector<double>& s)
        : origin(path.front().pos), arc(3 * (s.size() - 1) + 1),
          points(arc.size()), integrals(arc.size()) {
        const size_t last = s.size() - 1;
        const V3 end = path.back().pos - origin;
        for (size_t k = 0; k <= last; ++k) {
            const V3 p = path[k].pos - origin;
            arc[last + k] = s[k]; points[last + k] = p;
            // Odd reflection pins both process endpoints and continues their tangents.
            arc[last - k] = -s[k]; points[last - k] = p * -1;
            arc[3 * last - k] = 2 * s.back() - s[k];
            points[3 * last - k] = end * 2 - p;
        }
        for (size_t k = 1; k < arc.size(); ++k) integrals[k] = integrate(k - 1, arc[k]);
    }

    V3 average(double s, double radius) const {
        const double lo = s - radius, hi = s + radius;
        const auto a = integralAt(lo), b = integralAt(s), c = integralAt(hi);
        return origin + ((b.moment - a.moment) - (b.area - a.area) * lo
            + (c.area - b.area) * hi - (c.moment - b.moment)) * (1 / (radius * radius));
    }

private:
    struct Integral { V3 area, moment; };
    Integral integralAt(double s) const {
        const auto pos = std::upper_bound(arc.begin(), arc.end(), s) - arc.begin();
        return integrate(pos > 0 ? std::min(arc.size() - 2, static_cast<size_t>(pos - 1)) : 0, s);
    }
    Integral integrate(size_t i, double s) const {
        const double t = arc[i], h = s - t;
        const V3 p = points[i], v = (points[i + 1] - p) * (1 / (arc[i + 1] - t));
        const V3 area = p * h + v * (h * h / 2);
        return {integrals[i].area + area,
                integrals[i].moment + area * t + p * (h * h / 2) + v * (h * h * h / 3)};
    }
    V3 origin;
    std::vector<double> arc;
    std::vector<V3> points;
    std::vector<Integral> integrals;
};

double maxPathTurn(const std::vector<CornerPoint>& run) {
    double largest = 0;
    for (size_t k = 1; k + 1 < run.size(); ++k) {
        const V3 a = run[k].pos - run[k - 1].pos;
        const V3 b = run[k + 1].pos - run[k].pos;
        if (norm(a) > 1e-6 && norm(b) > 1e-6)
            largest = std::max(largest, std::atan2(norm(cross(a, b)), dot(a, b)) * kRad);
    }
    return largest;
}

template <typename Accept>
bool roundAndSampleTrajectory(std::vector<CornerPoint>& points, double internalStep, double outputStep,
                              double window, double standoff, Accept accept) {
    std::vector<CornerPoint> assembled, dense, candidate;
    assembled.reserve(points.size());
    int rounded = 0, restored = 0;
    bool changed = false;
    auto appendRun = [&](auto&& append, size_t begin, size_t end) -> void {
        std::vector<double> arc{0.0};
        for (size_t i = begin + 1; i < end; ++i) {
            const double gap = norm(points[i].pos - points[i - 1].pos);
            arc.push_back(arc.back() + gap);
        }
        const std::vector<CornerPoint> run(points.begin() + begin, points.begin() + end);
        if (run.size() < 2) {
            assembled.push_back(run.front());
            return;
        }
        const double oldTurn = maxPathTurn(run);
        const IntegratedPath curve(run, arc);
        const double radius = std::min(window, 0.5 * arc.back());
        const bool canSplit = end - begin >= 12;
        int failedAt = -1;
        for (double scale : {1.0, 0.75, 0.5, 0.25, 0.0}) {
            if (scale == 0.0 && canSplit && oldTurn > 35.0) break;
            const double reach = radius * scale;
            dense.clear();
            if (scale == 0.0 || run.size() < 3) {
                dense = run;
            } else {
                const double ds = std::clamp(std::min(outputStep, reach) / 4.0, 0.125, 1.0);
                const size_t steps = static_cast<size_t>(std::ceil(arc.back() / ds));
                dense.push_back(run.front());
                auto sourceAt = [&](double s) {
                    const auto index = std::upper_bound(arc.begin(), arc.end(), s) - arc.begin();
                    const size_t j = std::min(arc.size() - 2, static_cast<size_t>(index - 1));
                    const double t = (s - arc[j]) / (arc[j + 1] - arc[j]);
                    const auto& a = run[j];
                    const auto& b = run[j + 1];
                    return CornerPoint{a.pos + (b.pos - a.pos) * t, blendAxis(a.gun, b.gun, t),
                        a.speed, a.spray, a.kind, t < 0.5 ? a.source : b.source};
                };
                for (size_t k = 1; k < steps; ++k) {
                    const double s = arc.back() * k / steps;
                    const auto raw = sourceAt(s);
                    const V3 p = curve.average(s, reach);
                    if (norm(p - raw.pos) > 0.08 * standoff) { dense.clear(); break; }
                    double fade = std::clamp(std::min(s, arc.back() - s) / reach, 0.0, 1.0);
                    fade = fade * fade * (3 - 2 * fade);
                    const V3 mean = blendAxis(sourceAt(std::max(0.0, s - reach * 0.5)).gun,
                        sourceAt(std::min(arc.back(), s + reach * 0.5)).gun, 0.5);
                    const double turn = std::atan2(norm(cross(raw.gun, mean)), dot(raw.gun, mean)) * kRad;
                    const V3 axis = blendAxis(raw.gun, mean, turn > 3 ? std::min(fade, 3 / turn) : fade);
                    // Maintain the nominal spray footprint while moving the TCP into the corner.
                    const V3 gun = raw.spray ? unit(raw.pos + axis * standoff - p) : axis;
                    dense.push_back({p, gun, raw.speed, raw.spray, raw.kind, raw.source});
                }
                if (dense.empty()) continue;
                dense.push_back(run.back());
            }
            if (!sampleAlongCurve(dense, outputStep, candidate)) continue;
            if (maxPathTurn(candidate) > oldTurn + 1.0) continue;
            const int failure = accept(candidate);
            if (failure >= 0) {
                failedAt = failure;
                continue;
            }
            if (scale > 0 && run.size() >= 3) ++rounded;
            changed |= scale > 0 || candidate.size() != run.size();
            assembled.insert(assembled.end(), candidate.begin(), candidate.end());
            return;
        }
        if (canSplit) {
            // Restrict a failed geometric check to its neighborhood. For a sampling
            // failure, split near the middle where the path is closest to straight.
            const size_t guard = std::max<size_t>((end - begin) / 4, 4);
            const size_t center = failedAt >= 0
                ? std::clamp(static_cast<size_t>(failedAt), begin + guard, end - guard)
                : begin + (end - begin) / 2;
            const size_t half = std::min<size_t>((end - begin) / 8, 8);
            const size_t first = center - std::min(half, center - begin - 4);
            const size_t last = std::min(end - 4, center + half);
            size_t split = center;
            double best = 1e300;
            for (size_t i = first; i <= last; ++i) {
                const V3 edge = unit(points[i].pos - points[i - 1].pos);
                const V3 a = unit(points[i - 1].pos - points[i - 2].pos);
                const V3 b = unit(points[i + 1].pos - points[i].pos);
                const double score = 2.0 - dot(a, edge) - dot(edge, b);
                if (score + 1e-6 < best) { split = i; best = score; }
            }
            append(append, begin, split);
            append(append, split, end);
            return;
        }
        ++restored;
        assembled.insert(assembled.end(), run.begin(), run.end());
    };
    for (size_t begin = 0; begin < points.size();) {
        size_t end = begin + 1;
        for (; end < points.size(); ++end) {
            const double gap = norm(points[end].pos - points[end - 1].pos);
            if (!sameProcess(points[end - 1], points[end]) || gap < 1e-6 ||
                gap > std::max(20.0, 4.0 * internalStep)) break;
        }
        appendRun(appendRun, begin, end);
        begin = end;
    }
    std::fprintf(stderr, "[转角平滑] %d 段已圆滑，%d 段保留原路径 · 最大方向变化 %.2f° → %.2f° · 点距 %.2f mm · %zu → %zu 点\n",
                 rounded, restored, maxPathTurn(points), maxPathTurn(assembled),
                 outputStep, points.size(), assembled.size());
    if (maxPathTurn(assembled) > maxPathTurn(points) + 1.0) return false;
    points.swap(assembled);
    return changed;
}

} // namespace
