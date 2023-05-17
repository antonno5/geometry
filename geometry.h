#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

struct Point {
    double x = 0, y = 0;
    static constexpr const double EPS = 1e-10;

    static bool eps_equal(double first, double second) {
        return std::abs(first - second) < EPS;
    }

    static bool eps_less(double first, double second) {
        return first < second + EPS;
    }

    static double to_radians(double angle) {
        const double grad = 180;
        return angle / grad * acos(-1);
    }

    Point() = default;

    Point(double in_x, double in_y) : x(in_x), y(in_y) {}

    double length() const {
        return sqrt(x * x + y * y);
    }

    bool operator==(Point that) const {
        return eps_equal(x, that.x) && eps_equal(y, that.y);
    }

    bool operator!=(Point that) const {
        return !(*this == that);
    }

    Point& operator+=(Point that) {
        x += that.x;
        y += that.y;
        return *this;
    }

    friend Point operator+(Point first, Point second) {
        return first += second;
    }

    Point operator-() const {
        return {-x, -y};
    }

    Point& operator-=(Point that) {
        return *this += -that;
    }

    friend Point operator-(Point first, Point second) {
        return first -= second;
    }

    Point& operator*=(double coefficient) {
        x *= coefficient;
        y *= coefficient;
        return *this;
    }

    friend Point operator*(Point point, double coefficient) {
        return point *= coefficient;
    }

    Point& operator/=(double coefficient) {
        x /= coefficient;
        y /= coefficient;
        return *this;
    }

    friend Point operator/(Point point, double coefficient) {
        return point /= coefficient;
    }

    double operator*(Point that) const {
        return x * that.x + y * that.y;
    }

    double operator%(Point that) const {
        return x * that.y - that.x * y;
    }

    friend bool collinear(Point first, Point second) {
        return eps_equal(first % second, 0);
    }

    Point rotate(double angle) const {
        double x1 = x, y1 = y;
        return {cos(angle) * x1 - sin(angle) * y1, sin(angle) * x1 + cos(angle) * y1};
    }

    Point rotate(Point center, double angle) const {
        return center + (*this - center).rotate(angle);
    }

    Point normalized() const {
        return {x / this->length(), y / this->length()};
    }

    Point reflect(Point center) const {
        return center - (*this - center);
    }

    Point scale(Point center, double coefficient) const {
        return center + (*this - center) * coefficient;
    }

    static Point relation(Point first, Point second, double coefficient) {
        return first + (second - first) * coefficient;
    }

    static Point middle(Point first, Point second) {
        return (first + second) / 2;
    }

    static double get_angle(Point first, Point second) {
        double cos_angle = first * second / first.length() / second.length();
        if (cos_angle < -1)
            cos_angle = -1;
        if (cos_angle > 1)
            cos_angle = 1;
        return acos(cos_angle) * (second % first > 0 ? 1 : -1);
    }

    static Point norm(Point first, Point second) {
        return (second - first).rotate(-acos(-1) / 2).normalized();
    }

    bool inSegment(Point first, Point second) {
        return eps_less((*this - first).length(), (second - first).length()) &&
                eps_less((*this - second).length(), (second - first).length());
    }
};

class Line {
    Point m_point, m_dir;
public:
    Point point() {
        return m_point;
    }

    Point dir() {
        return m_dir;
    }

    Line() = default;

    Line(Point point1, Point point2) : m_point(point1), m_dir(point2 - point1) {}

    Line(double coefficient, double shift) : m_point(0, shift), m_dir(1, coefficient) {}

    Line(Point point, double coefficient) : m_point(point),
                                      m_dir(1, coefficient) {}

    static Line construct(Point point, Point dir) {
        Line answer;
        answer.m_point = point;
        answer.m_dir = dir;
        return answer;
    }

    Line shift(Point new_point) const {
        Line answer = *this;
        answer.m_point = new_point;
        return answer;
    }

    bool operator==(const Line& that) const {
        return collinear(m_dir, that.m_dir) && collinear(m_point - that.m_point, m_dir);
    }

    bool operator!=(const Line& that) const {
        return !(*this == that);
    }

    Line rotate(double angle) const {
        return Line::construct(m_point, m_dir.rotate(angle));
    }

    static Point intersection(const Line& first, const Line& second) {
        double coeff1 = second.m_point % second.m_dir - first.m_point % second.m_dir;
        double coeff2 = first.m_dir % second.m_dir;
        return first.m_point + first.m_dir * (coeff1 / coeff2);
    }

    static Point perpendicular(const Line& line, Point point) {
        return intersection(line, line.rotate(acos(-1) / 2).shift(point));
    }

    static Point reflect(const Line& line, Point point) {
        return point + (perpendicular(line, point) - point) * 2;
    }

    static Line midperpendicular(Point first, Point second) {
        return construct(Point::middle(first, second), Point::norm(first, second));
    }
};

class Shape {
public:
    virtual bool operator==(const Shape& that) const = 0;

    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool isCongruentTo(const Shape& that) const = 0;

    virtual bool isSimilarTo(const Shape& that) const = 0;

    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;

    virtual void reflect(const Point& center) = 0;

    virtual void reflect(const Line& axis) = 0;

    virtual void scale(const Point& center, double coefficient) = 0;

    virtual bool operator!=(const Shape& that) const = 0;

    virtual ~Shape() = default;
};

class Polygon : public Shape{
protected:
    std::vector<Point> vertices;

    Point at(int i) const {
        return vertices[static_cast<size_t>(i % static_cast<int>(vertices.size()) + static_cast<int>(vertices.size()))
        % vertices.size()];
    }

    Point at(size_t i) const {
        return vertices[i % vertices.size()];
    }

    bool isSimilarInOneWay(const Polygon& that, int i, bool invert) const {
        bool equal1 = true, equal2 = true;
        int inv = (invert ? -1 : 1);
        double coefficient = (at(1) - at(0)).length() / (that.at(i + inv) - that.at(i)).length();
        for (int j = 0; j < static_cast<int>(verticesCount()); ++j) {
            auto this_first = at(j + 1) - at(j);
            auto this_second = at(j + 2) - at(j + 1);
            auto that_first = that.at(i + (j + 1) * inv) - that.at(i + j * inv);
            auto that_second = that.at(i + (j + 2) * inv) - that.at(i + (j + 1) * inv);
            if (!Point::eps_equal(this_first.length() / that_first.length(), coefficient)) {
                equal1 = equal2 = false;
            }
            if (!Point::eps_equal(Point::get_angle(this_second, this_first),
                                  Point::get_angle(that_second, that_first))) {
                equal1 = false;
            }
            if (!Point::eps_equal(Point::get_angle(this_second, this_first),
                                  -Point::get_angle(that_second, that_first))) {
                equal2 = false;
            }
        }
        return equal1 || equal2;
    }

    bool equalInOneWay(const Polygon& that, int ind, bool invert) const {
        bool equal = true;
        int inv = (invert ? -1 : 1);
        for (int i = 0; i < static_cast<int>(that.verticesCount()); ++i) {
                if (at(i) != that.at(ind + i * inv))
                    equal = false;
        }
        return equal;
    }

    double area_signed() const {
        double ans = (vertices.front() % vertices.back());
        for (size_t i = 0; i < vertices.size() - 1; ++i)
            ans += vertices[i + 1] % vertices[i];
        return ans;
    }
public:
    size_t verticesCount() const {
        return vertices.size();
    }

    std::vector<Point> getVertices() const {
        return vertices;
    }

    bool isConvex() const {
        bool equal1 = true, equal2 = true;
        for (size_t i = 0; i < vertices.size(); ++i) {
            const Point& first = at(i);
            const Point& second = at(i + 1);
            const Point& third = at(i + 2);
            if (Point::eps_less((third - second) % (second - first), 0)) {
                equal1 = false;
            }
            if (Point::eps_less(0, (third - second) % (second - first))) {
                equal2 = false;
            }
        }
        return equal1 || equal2;
    }

    Polygon() = default;

    double area() const override {
        return std::abs(area_signed()) / 2;
    }

    Polygon(std::vector<Point> in_vertices) : vertices(in_vertices) {}

    template<typename ...T>
    explicit Polygon(T... points) {
        ((vertices.push_back(points)), ...);
    }

    explicit Polygon(std::initializer_list<Point> list) {
        for (auto i : list) {
            vertices.push_back(i);
        }
    }

    double perimeter() const override {
        double ans = (vertices.front() - vertices.back()).length();
        for (size_t i = 0; i < vertices.size() - 1; ++i)
            ans += (vertices[i + 1] - vertices[i]).length();
        return ans;
    }

    bool operator==(const Polygon& that) const {
        if (verticesCount() != that.verticesCount())
            return false;
        size_t ind = that.verticesCount();
        for (size_t i = 0; i < that.verticesCount(); ++i) {
            if (vertices.front() == that.vertices[i])
                ind = i;
        }
        return (ind != that.verticesCount()) && (equalInOneWay(that, static_cast<int>(ind), false)
        || equalInOneWay(that, static_cast<int>(ind), true));
    }

    bool operator==(const Shape& that) const override {
        auto polygon_that_pointer = dynamic_cast<const Polygon*>(&that);
        if (polygon_that_pointer == nullptr)
            return false;
        return *this == *polygon_that_pointer;
    }

    bool operator!=(const Polygon& that) const {
        return !(*this == that);
    }

    bool operator!=(const Shape& that) const override {
        return !(*this == that);
    }

    bool isSimilarTo(const Shape& shape_that) const override {
        auto *that_pointer = dynamic_cast<const Polygon*>(&shape_that);
        if (that_pointer == nullptr)
            return false;
        if (verticesCount() != that_pointer->verticesCount())
            return false;
        for (size_t i = 0; i < verticesCount(); ++i) {
            if (isSimilarInOneWay(*that_pointer, static_cast<int>(i), false)
                || isSimilarInOneWay(*that_pointer, static_cast<int>(i), true))
                return true;
        }
        return false;
    }

    bool isCongruentTo(const Shape& shape_that) const override {
        auto *that_pointer = dynamic_cast<const Polygon*>(&shape_that);
        if (that_pointer == nullptr)
            return false;
        return this->isSimilarTo(*that_pointer) && Point::eps_equal(area(), that_pointer->area());
    }

    bool containsPoint(const Point& point) const override {
        double min_len = -1;
        int ind = -1;
        bool is_vertice = true;
        Point vector;
        for (size_t i = 0; i < vertices.size(); ++i) {
            auto per = at(i) - point;
            if (per.length() < min_len || ind == -1) {
                min_len = per.length();
                vector = per;
                ind = static_cast<int>(i);
            }
        }
        for (size_t i = 0; i < vertices.size(); ++i) {
            auto first = at(i);
            auto second = at(i + 1);
            auto per_point = Line::perpendicular(Line(first, second), point);
            auto per = per_point - point;
            if (per_point.inSegment(first, second) && per.length() < min_len) {
                min_len = per.length();
                ind = static_cast<int>(i);
                vector = per;
                is_vertice = false;
            }
        }
        Point vector_norm;
        if (is_vertice) {
            auto first = at(ind - 1);
            auto second = at(ind);
            auto third = at(ind + 1);
            vector_norm = -((first - second).normalized() + (third - second).normalized()) / 2;
            if (Point::eps_less((third - second) % (second - first), 0)) {
                vector_norm = -vector_norm;
            }
        } else {
            auto first = at(ind);
            auto second = at(ind + 1);
            vector_norm = -Point::norm(first, second);
        }
        if (area_signed() < 0)
            vector_norm = -vector_norm;
        return Point::eps_less(0, vector * vector_norm);
    }

    void rotate(const Point& center, double angle) override {
        angle = Point::to_radians(angle);
        for (auto & vertice : vertices)
            vertice = vertice.rotate(center, angle);
    }

    void reflect(const Point& center) override {
        for (auto & vertice : vertices)
            vertice = vertice.reflect(center);
    }

    void reflect(const Line& axis) override {
        for (auto & vertice : vertices)
            vertice = Line::reflect(axis, vertice);
    }

    void scale(const Point& center, double coefficient) override {
        for (auto & vertice : vertices)
            vertice = vertice.scale(center, coefficient);
    }
};

class Ellipse : public Shape{
protected:
    Point focus1, focus2;
    double dist;
public:
    Ellipse(Point in_focus1, Point in_focus2, double in_dist) : focus1(in_focus1), focus2(in_focus2), dist(in_dist) {}

    std::pair<Point, Point> focuses() const {
        return {focus1, focus2};
    }

    Point center() const {
        return Point::middle(focus1, focus2);
    }

    double get_a() const {
        return dist / 2;
    }

    double get_b() const {
        return sqrt(get_a() * get_a() - (focus1 - center()).length() * (focus1 - center()).length());
    }

    double eccentricity() const {
        return sqrt(1 - (get_b() * get_b()) / (get_a() * get_a()));
    }

    std::pair<Line, Line> directrices() const {
        Point dir = Point::norm(focus1, center());
        Point vector = (focus1 - center()).normalized() * get_a() / eccentricity();
        return {Line::construct(center() + vector, dir), Line::construct(center() - vector, dir)};
    }

    double perimeter() const override {
        const double number = 3;
        return acos(-1) * (number * (get_a() + get_b()) - sqrt((number * get_a() + get_b()) *
        (get_a() + number * get_b())));
    }

    double area() const override {
        return acos(-1) * get_a() * get_b();
    }

    bool operator==(const Ellipse& that) const {
        return Point::eps_equal(dist, that.dist) && (((focus1 == that.focus1) && (focus2 == that.focus2))
                                                           || ((focus1 == that.focus2) && (focus2 == that.focus1)));
    }

    bool operator==(const Shape& that) const override {
        auto ellipse_that_pointer = dynamic_cast<const Ellipse*>(&that);
        if (ellipse_that_pointer == nullptr)
            return false;
        return *this == *ellipse_that_pointer;
    }

    bool operator!=(const Ellipse& that) const {
        return !(*this == that);
    }

    bool operator!=(const Shape& that) const override {
        return !(*this == that);
    }

    bool isCongruentTo(const Shape& shape_that) const override {
        auto *that_pointer = dynamic_cast<const Ellipse*>(&shape_that);
        if (that_pointer == nullptr)
            return false;
        return Point::eps_equal(dist, that_pointer->dist) && Point::eps_equal((focus1 - focus2).length(),
                                                                     (that_pointer->focus1 - that_pointer->focus2).length());
    }

    bool isSimilarTo(const Shape& shape_that) const override {
        auto *that_pointer = dynamic_cast<const Ellipse*>(&shape_that);
        if (that_pointer == nullptr)
            return false;
        return Point::eps_equal(eccentricity(), that_pointer->eccentricity());
    }

    bool containsPoint(const Point& point) const override {
        return Point::eps_less((point - focus1).length() + (point - focus2).length(), dist);
    }

    void rotate(const Point& center, double angle) override {
        angle = Point::to_radians(angle);
        focus1.rotate(center, angle);
        focus2.rotate(center, angle);
    }

    void reflect(const Point& center) override {
        focus1 = focus1.reflect(center);
        focus2 = focus2.reflect(center);
    }

    void reflect(const Line& axis) override {
        focus1 = Line::reflect(axis, focus1);
        focus2 = Line::reflect(axis, focus2);
    }

    void scale(const Point& center, double coefficient) override {
        focus1 = focus1.scale(center, coefficient);
        focus2 = focus2.scale(center, coefficient);
        dist *= coefficient;
    }
};

class Circle : public Ellipse{
public:
    double radius() {
        return dist;
    }

    Circle(Point center, double radius) : Ellipse(center, center, radius) {}
};

class Rectangle : public Polygon{
public:
    Point center() {
        return Point::middle(vertices[0], vertices[2]);
    }

    std::pair<Line, Line> diagonals() {
        const size_t magic_number = 3;
        return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[magic_number])};
    }

    Rectangle(Point first, Point second, double coefficient) {
        if (coefficient < 1)
            coefficient = 1 / coefficient;
        double angle = atan(coefficient);
        Point rotated = (Point::middle(first, second) - first).rotate(angle);
        double length = (first - second).length() / sqrt(coefficient * coefficient + 1);
        rotated = rotated.normalized() * length;
        vertices = {first, first + rotated, second, second - rotated};
    }

    Rectangle(Point first, Point second, Point third, Point fourth) : Polygon(first, second, third, fourth) {}
};

class Square : public Rectangle{
public:
    Circle circumscribedCircle() {
        return {Point::middle(vertices[0], vertices[2]), (vertices[0] - center()).length()};
    }

    Circle inscribedCircle() {
        return {Point::middle(vertices[0], vertices[2]), (vertices[0] - vertices[1]).length() / 2};
    }

    Square(Point first, Point second) : Rectangle(first, second, 1) {}

    Square(Point first, Point second, Point third, Point fourth) : Rectangle(first, second, third, fourth) {}
};

class Triangle : public Polygon{
public:
    Triangle(Point first, Point second, Point third) : Polygon(first, second, third) {}

    Circle circumscribedCircle() const {
        Line MidPerpendicular1 = Line::midperpendicular(vertices[0], vertices[1]);
        Line MidPerpendicular2 = Line::midperpendicular(vertices[1], vertices[2]);
        Point center = Line::intersection(MidPerpendicular1, MidPerpendicular2);
        return {center, (vertices[0] - center).length()};
    }

    Circle inscribedCircle() const {
        double relation1 = (vertices[1] - vertices[0]).length() /
                           ((vertices[1] - vertices[0]).length() + (vertices[2] - vertices[0]).length());
        Line Bisect1 = Line::construct(vertices[0],
                                       Point::relation(vertices[1], vertices[2], relation1) - vertices[0]);
        double relation2 = (vertices[0] - vertices[1]).length() /
                           ((vertices[0] - vertices[1]).length() + (vertices[2] - vertices[1]).length());
        Line Bisect2 = Line::construct(vertices[1],
                                       Point::relation(vertices[0], vertices[2], relation2) - vertices[1]);
        Point center = Line::intersection(Bisect1, Bisect2);
        return {center, area() / perimeter() * 2};
    }

    Point centroid() const {
        return (vertices[0] + vertices[1] + vertices[2]) / static_cast<double>(vertices.size());
    }

    Point orthocenter() const {
        Line height1(vertices[0], Line::perpendicular(Line(vertices[1], vertices[2]), vertices[0]));
        Line height2(vertices[1], Line::perpendicular(Line(vertices[0], vertices[2]), vertices[1]));
        return Line::intersection(height1, height2);
    }

    Circle ninePointsCircle() const {
        return Triangle(Point::middle(vertices[0], vertices[1]), Point::middle(vertices[1], vertices[2]),
                        Point::middle(vertices[2], vertices[0])).circumscribedCircle();
    }

    Line EulerLine() const {
        return {ninePointsCircle().center(), centroid()};
    }
};

