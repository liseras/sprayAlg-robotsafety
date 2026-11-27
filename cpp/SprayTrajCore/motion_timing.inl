// Position and gun axis stay on the safety layer's output segments. Time laws
// for TCP travel and each measured joint share the same segment durations.
namespace {

struct MotionSpan {
    double length = 0, gunTurn = 0, cap = 0, duration = 0;
    bool jointChecked = false;
};

template <size_t N>
std::array<double, N - 1> bezierDerivative(const std::array<double, N>& curve, double duration) {
    std::array<double, N - 1> derivative{};
    for (size_t i = 0; i + 1 < N; ++i)
        derivative[i] = (curve[i + 1] - curve[i]) * (N - 1) / duration;
    return derivative;
}

// Two subdivisions tighten the convex-hull bound while keeping it conservative.
template <size_t N>
double bezierMagnitudeBound(std::array<double, N> curve, int subdivisions = 2) {
    if (subdivisions == 0) {
        double bound = 0;
        for (double v : curve) bound = std::max(bound, std::fabs(v));
        return bound;
    }
    std::array<double, N> left{}, right{};
    left[0] = curve[0];
    right[N - 1] = curve[N - 1];
    for (size_t j = 1; j < N; ++j) {
        for (size_t i = 0; i < N - j; ++i) curve[i] = (curve[i] + curve[i + 1]) * 0.5;
        left[j] = curve[0];
        right[N - 1 - j] = curve[N - 1 - j];
    }
    return std::max(bezierMagnitudeBound(left, subdivisions - 1),
                    bezierMagnitudeBound(right, subdivisions - 1));
}

struct MotionPeaks { double speed, accel, jerk; };

MotionPeaks motionPeaks(double delta, double firstSpeed, double lastSpeed, double duration,
                       double firstAcceleration = 0, double lastAcceleration = 0) {
    const double a = firstSpeed * duration / 5.0, b = lastSpeed * duration / 5.0;
    const double scale = duration * duration / 20;
    const std::array<double, 6> curve{0, a, 2 * a + firstAcceleration * scale,
        delta - 2 * b + lastAcceleration * scale, delta - b, delta};
    const auto speed = bezierDerivative(curve, duration);
    const auto accel = bezierDerivative(speed, duration);
    const auto jerk = bezierDerivative(accel, duration);
    return {bezierMagnitudeBound(speed), bezierMagnitudeBound(accel), bezierMagnitudeBound(jerk)};
}

// Keep the acceleration-dependent control points ordered. Split the available
// middle control interval between both ends so they can be limited independently.
std::pair<double, double> jointAccelerationRange(double delta, double startSpeed, double endSpeed,
                                                  double duration, bool atEnd) {
    const double room = std::max(0.0, std::fabs(delta) -
        0.4 * duration * (std::fabs(startSpeed) + std::fabs(endSpeed)));
    const double middle = 10 * room / (duration * duration);
    const double near = 4 * std::fabs(atEnd ? endSpeed : startSpeed) / duration;
    const double lo = atEnd ? -middle : -near;
    const double hi = atEnd ? near : middle;
    return delta >= 0 ? std::make_pair(lo, hi) : std::make_pair(-hi, -lo);
}

} // namespace

int robot_safety_plan_motion(const RobotSafetyPoint* points, int count,
                             const RobotMotionLimits* limits, RobotMotionPoint* timed) {
    if (!points || !limits || !timed || count < 2) return 0;
    const RobotMotionLimits& l = *limits;
    for (double value : {l.straightSpeedMmSec, l.turnSpeedMmSec, l.tcpAccelMmSec2,
                         l.tcpJerkMmSec3, l.gunAxisSpeedDegSec, l.jointSpeedDegSec,
                         l.jointAccelDegSec2, l.jointJerkDegSec3}) {
        if (!std::isfinite(value) || value <= 0) return 0;
    }

    constexpr double stationary = 1e-6;
    constexpr unsigned invalidJoint = RS_VIOL_IK_FAIL | RS_VIOL_LINK_HIT | RS_VIOL_TOOL_HIT |
                                      RS_VIOL_SELF_COLLISION;
    const double turnSpeed = std::min(l.turnSpeedMmSec, l.straightSpeedMmSec);
    std::vector<MotionSpan> spans(count);     // Each entry describes travel into its point.
    std::vector<double> reference(count), nodeCap(count), cornerCap(count, l.straightSpeedMmSec);
    std::vector<bool> stop(count, false);
    stop.front() = stop.back() = true;
    bool moving = false;
    auto pos = [&](int i) { return V3(points[i].x, points[i].y, points[i].z); };
    auto axis = [&](int i) { return V3(points[i].i, points[i].j, points[i].k); };
    auto validJoints = [invalidJoint](const RobotSafetyPoint& point) {
        if (!point.ikOk || point.marginDeg < 0 || (point.violation & invalidJoint)) return false;
        return std::all_of(point.theta, point.theta + 6, [](double q) { return std::isfinite(q); });
    };
    for (int i = 0; i < count; ++i) {
        timed[i] = {};
        const V3 p = pos(i), a = axis(i);
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
            !std::isfinite(norm(a)) || norm(a) <= stationary ||
            !std::isfinite(points[i].speed)) return 0;
    }
    for (int i = 1; i < count; ++i) {
        MotionSpan& s = spans[i];
        s.length = norm(pos(i) - pos(i - 1));
        s.gunTurn = angleBetweenDeg(axis(i - 1), axis(i));
        const auto& a = points[i - 1];
        const auto& b = points[i];
        double jointDelta = 0;
        if (validJoints(a) && validJoints(b)) {
            for (int j = 0; j < 6; ++j)
                jointDelta = std::max(jointDelta, std::fabs(b.theta[j] - a.theta[j]));
            s.jointChecked = jointDelta <= 30.0;
            if (!s.jointChecked) {
                timed[i].flags |= RS_MOTION_JOINT_BRANCH_JUMP;
                stop[i - 1] = stop[i] = true;
            }
        }
        if (!s.jointChecked) timed[i].flags |= RS_MOTION_JOINT_UNCHECKED;
        if (s.length <= stationary) stop[i - 1] = stop[i] = true;
        moving |= s.length > stationary || s.gunTurn > 1e-5 || (s.jointChecked && jointDelta > 1e-5);
        if (s.gunTurn > 12.0 || (timed[i].flags & RS_MOTION_JOINT_BRANCH_JUMP)) {
            cornerCap[i - 1] = std::min(cornerCap[i - 1], turnSpeed);
            cornerCap[i] = std::min(cornerCap[i], turnSpeed);
        }
        const double nominalA = a.speed > 0 ? a.speed : l.straightSpeedMmSec;
        const double nominalB = b.speed > 0 ? b.speed : l.straightSpeedMmSec;
        s.cap = std::min({nominalA, nominalB, l.straightSpeedMmSec});
    }
    if (!moving) return 0;
    for (int i = 1; i + 1 < count; ++i) {
        if (spans[i].length <= stationary || spans[i + 1].length <= stationary) continue;
        const double turn = angleBetweenDeg(pos(i) - pos(i - 1), pos(i + 1) - pos(i));
        if (turn >= 45.0) {
            stop[i] = true;
            timed[i].flags |= RS_MOTION_STOP_TURN;
        }
        const double weight = clampd((turn - 5.0) / 15.0, 0.0, 1.0);
        const double cap = l.straightSpeedMmSec * (1 - weight) + turnSpeed * weight;
        cornerCap[i] = std::min(cornerCap[i], cap);
        if (turn > 1e-6) {
            const double radius = std::min(spans[i].length, spans[i + 1].length) /
                                  (2 * std::sin(turn * kDeg * 0.5));
            cornerCap[i] = std::min(cornerCap[i], std::sqrt(l.tcpAccelMmSec2 * radius));
        }
    }
    for (int i = 1; i < count; ++i) {
        MotionSpan& s = spans[i];
        s.cap = std::min({s.cap, cornerCap[i - 1], cornerCap[i]});
        if (s.length > stationary && s.gunTurn > 1e-6)
            s.cap = std::min(s.cap, s.length * l.gunAxisSpeedDegSec / s.gunTurn);
        timed[i].incomingSpeedCapMmSec = s.cap;
    }
    for (int i = 1; i + 1 < count; ++i) {
        if (!stop[i]) nodeCap[i] = std::min(spans[i].cap, spans[i + 1].cap);
    }

    // With duration 2*ds/(v0+v1), the quintic devolves to a cubic smoothstep
    // in speed. Its peak acceleration is 1.5*dv/dt and peak jerk is 6*dv/dt².
    // Use the larger endpoint speed to bound jerk monotonically in each pass.
    auto reachableSpeed = [&](double from, double target, double distance) {
        if (target <= from) return target;
        const double high = std::min(target, std::sqrt(from * from + l.tcpAccelMmSec2 * distance * 1.333));
        const double budget = l.tcpJerkMmSec3 * distance * distance / 6 * 0.999;
        if ((high - from) * high * high <= budget) return high;
        double lo = from, hi = high;
        for (int j = 0; j < 24; ++j) {
            const double v = (lo + hi) * 0.5;
            if ((v - from) * v * v <= budget) lo = v;
            else hi = v;
        }
        return lo;
    };
    const MotionPeaks fromRest = motionPeaks(1, 0, 0, 1);
    auto buildProfile = [&] {
        reference = nodeCap;
        for (int i = 1; i < count; ++i)
            reference[i] = reachableSpeed(reference[i - 1], reference[i], spans[i].length);
        for (int i = count - 1; i > 0; --i)
            reference[i - 1] = reachableSpeed(reference[i], reference[i - 1], spans[i].length);
        for (int i = 1; i < count; ++i) {
            MotionSpan& s = spans[i];
            const double sum = reference[i - 1] + reference[i];
            if (sum > stationary) {
                s.duration = 2 * s.length / sum;
                continue;
            }
            // Both end speeds are zero: independently time the small move or
            // in-place rotation, using the rest-to-rest quintic.
            double jointDelta = 0;
            if (s.jointChecked) {
                for (int j = 0; j < 6; ++j)
                    jointDelta = std::max(jointDelta, std::fabs(points[i].theta[j] - points[i - 1].theta[j]));
            }
            s.duration = 1.001 * std::max({0.001, fromRest.speed * s.length / s.cap,
                std::sqrt(fromRest.accel * s.length / l.tcpAccelMmSec2),
                std::cbrt(fromRest.jerk * s.length / l.tcpJerkMmSec3),
                fromRest.speed * s.gunTurn / l.gunAxisSpeedDegSec,
                fromRest.speed * jointDelta / l.jointSpeedDegSec,
                std::sqrt(fromRest.accel * jointDelta / l.jointAccelDegSec2),
                std::cbrt(fromRest.jerk * jointDelta / l.jointJerkDegSec3)});
        }
        for (int i = 1; i + 1 < count; ++i) {
            const auto& a = spans[i];
            const auto& b = spans[i + 1];
            timed[i].tcpSpeedMmSec = reference[i];
            if (!a.jointChecked || !b.jointChecked || reference[i] <= stationary) {
                std::fill_n(timed[i].jointVelocityDegSec, 6, 0.0);
                continue;
            }
            for (int j = 0; j < 6; ++j) {
                const double q = points[i].theta[j];
                const double before = q - points[i - 1].theta[j], after = points[i + 1].theta[j] - q;
                const double beforeRate = before / a.duration, afterRate = after / b.duration;
                double rate = 0;
                if (std::signbit(beforeRate) == std::signbit(afterRate)) {
                    const double estimate = (beforeRate * b.duration + afterRate * a.duration) /
                                            (a.duration + b.duration);
                    const double cap = 1.25 * std::min(std::fabs(beforeRate), std::fabs(afterRate));
                    rate = std::copysign(std::min(std::fabs(estimate), cap), estimate);
                }
                timed[i].jointVelocityDegSec[j] = rate;
            }
        }
        for (int i = 1; i + 1 < count; ++i) {
            const auto& a = spans[i];
            const auto& b = spans[i + 1];
            if (!a.jointChecked || !b.jointChecked || reference[i] <= stationary) {
                std::fill_n(timed[i].jointAccelerationDegSec2, 6, 0.0);
                continue;
            }
            for (int j = 0; j < 6; ++j) {
                const double q = points[i].theta[j];
                const double before = q - points[i - 1].theta[j], after = points[i + 1].theta[j] - q;
                const double velocity = timed[i].jointVelocityDegSec[j];
                const auto left = jointAccelerationRange(before, timed[i - 1].jointVelocityDegSec[j],
                                                         velocity, a.duration, true);
                const auto right = jointAccelerationRange(after, velocity, timed[i + 1].jointVelocityDegSec[j],
                                                          b.duration, false);
                const double estimate = 2 * (after / b.duration - before / a.duration) /
                                        (a.duration + b.duration);
                timed[i].jointAccelerationDegSec2[j] =
                    std::clamp(estimate, std::max(left.first, right.first), std::min(left.second, right.second));
            }
        }
    };
    auto requiredStretch = [&](int i) {
        const auto& s = spans[i];
        const auto& a = timed[i - 1];
        auto& b = timed[i];
        const MotionPeaks tcp = motionPeaks(s.length, a.tcpSpeedMmSec, b.tcpSpeedMmSec, s.duration);
        const double gunRate = s.length > stationary ? s.gunTurn * tcp.speed / s.length
            : motionPeaks(s.gunTurn, 0, 0, s.duration).speed;
        double stretch = std::max({1.0, tcp.speed / s.cap, std::sqrt(tcp.accel / l.tcpAccelMmSec2),
                                  std::cbrt(tcp.jerk / l.tcpJerkMmSec3), gunRate / l.gunAxisSpeedDegSec});
        if (!s.jointChecked) return stretch;
        MotionPeaks joint{};
        for (int j = 0; j < 6; ++j) {
            const MotionPeaks q = motionPeaks(points[i].theta[j] - points[i - 1].theta[j],
                                              a.jointVelocityDegSec[j], b.jointVelocityDegSec[j], s.duration,
                                              a.jointAccelerationDegSec2[j], b.jointAccelerationDegSec2[j]);
            joint.speed = std::max(joint.speed, q.speed);
            joint.accel = std::max(joint.accel, q.accel);
            joint.jerk = std::max(joint.jerk, q.jerk);
        }
        b.incomingJointSpeedDegSec = joint.speed;
        b.incomingJointAccelDegSec2 = joint.accel;
        b.incomingJointJerkDegSec3 = joint.jerk;
        return std::max({stretch, joint.speed / l.jointSpeedDegSec,
                         std::sqrt(joint.accel / l.jointAccelDegSec2),
                         std::cbrt(joint.jerk / l.jointJerkDegSec3)});
    };

    // Tighten velocity caps only near a span whose joint derivatives are too
    // high. Rebuild the forward/backward speed ramps across those caps.
    for (int pass = 0; pass < 24; ++pass) {
        buildProfile();
        bool adjusted = false;
        for (int i = 1; i < count; ++i) {
            const double factor = requiredStretch(i);
            if (factor <= 1.000001) continue;
            if (reference[i - 1] + reference[i] <= stationary) return 0;
            nodeCap[i - 1] = std::min(nodeCap[i - 1], reference[i - 1] / (factor * 1.01));
            nodeCap[i] = std::min(nodeCap[i], reference[i] / (factor * 1.01));
            adjusted = true;
        }
        if (!adjusted) break;
    }
    buildProfile();
    for (int i = 1; i < count; ++i) {
        timed[i].timeSec = timed[i - 1].timeSec + spans[i].duration;
        if (!std::isfinite(timed[i].timeSec) || timed[i].timeSec <= timed[i - 1].timeSec ||
            requiredStretch(i) > 1.000001) return 0;
    }
    return 1;
}
