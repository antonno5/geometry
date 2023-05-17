#pragma once

// By Dan Klishch

#include "geometry.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "tiny_test.hpp"

using testing::make_test;
using testing::PrettyTest;
using testing::TestGroup;

static double double_error = 1e-5;

// ===== Polygon =====
struct PolygonData {
  std::initializer_list<Point> points;
  bool is_convex = false;
  double perimeter = 0.0;
  double area = 0.0;

  bool do_raster = false;
  double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
  size_t samples_x = 0, samples_y = 0;
  std::string image = "";

  std::vector<Point> getPointsVector() const {
    return {points};
  }

  std::vector<Point> getReversedPoints() const {
    auto result = getPointsVector();
    std::reverse(result.begin(), result.end());
    return result;
  }
};

const std::vector<Point> manual_test_points = {
    Point{0.8, 0.8},
    Point{2.6, 2.8},
    Point{0.6, 2.2},
};

const PolygonData polygon_collection[] = {
  {
    .points = {
      {0.8, 0.8},
      {3.4, 1.8},
      {2.6, 2.8},
      {1.4, 1.8},
      {0.6, 2.2},
    },
    .is_convex = false,
    .perimeter = 7.936993191477734,
    .area = 2.52,

    .do_raster = true,
    .x0 = 0, .y0 = 0, .x1 = 4, .y1 = 3,
    .samples_x = 41, .samples_y = 31,
    .image = R"((


                          #
                         ##
                        ####
                       ######
                      ########
                    ###########
      #            ############
       ##         ##############
       ####      ################
       ######   ##################
       ############################
       #########################
       ######################
       ####################
        ################
        ##############
        ###########
        ########
        ######
        ###
        #








))",
  },

  {
    .points = {
      {0.1, 0.1},
      {3.4, 1.8},
      {2.6, 2.8},
      {1.4, 1.8},
    },
    .is_convex = true,
    .perimeter = 8.694910477225286,
    .area = 2.7,

    .do_raster = true,
    .x0 = 0, .y0 = 0, .x1 = 4, .y1 = 3,
    .samples_x = 41, .samples_y = 31,
    .image = R"((


                          #
                         ##
                        ####
                       ######
                      ########
                    ###########
                   ############
                  ##############
                 ################
                ##################
              #####################
              ###################
             ##################
            #################
           ################
           ##############
          #############
         ############
        ###########
        #########
       ########
      #######
     ######
     ####
    ###
   ##
  #
 #

))",
  },

  {
    .points = {
      {0.4, 1.4},
      {0.6, 1.4},
      {0.6, 1.8},
      {1, 1.8},
      {1, 1.4},
      {1.4, 1.4},
      {1.6, 2},
      {1.8, 1.6},
      {2.4, 1.4},
      {1.8, 1},
      {2.2, 0.8},
      {1.8, 0.4},
      {0.4, 0.4},
    },
    .is_convex = false,
    .perimeter = 7.646133935109306,
    .area = 2,

    .do_raster = true,
    .x0 = 0, .y0 = 0, .x1 = 2.6, .y1 = 2.2,
    .samples_x = 27, .samples_y = 23,
    .image = R"((


                #
                #
      #####     ##
      #####    ###
      #####    ####
      #####    #######
    #####################
    ###################
    ##################
    ################
    ###############
    #################
    ###################
    ##################
    #################
    ################
    ###############




))",
  },

  {
    .points = {
      {2, 1.8},
      {3.4, 1.8},
      {2.6, 2.8},
      {1, 1.4},
    },
    .is_convex = false,
    .perimeter = 5.8836869714604,
    .area = 1.08,
  },

  {
    .points = {
      {0, 0},
      {1, 0},
      {1, 1},
      {0, 1},
      {-0.5, 0.5},
    },
    .is_convex = true,
    .perimeter = 4.414213562373095,
    .area = 1.25,
  },

  {
    .points = {
      {0, 0},
      {1, 0},
      {1, 1},
      {0, 1},
      {0.5, 0.5},
    },
    .is_convex = false,
    .perimeter = 4.414213562373095,
    .area = 0.75,
  },

  {
    .points = {
      {0, 0},
      {1, 0},
      {0.5, 0.8660254037844386},
    },
    .is_convex = true,
    .perimeter = 3,
    .area = 0.4330127018922193,
  },

  {
    .points = {
      {1, 1},
      {1, 4},
      {5, 4},
    },
    .is_convex = true,
    .perimeter = 12,
    .area = 6,
  },

  {
    .points = {
      {2.2, 1.2},
      {3.4, 1.8},
      {1.2, 2},
    },
    .is_convex = true,
    .perimeter = 4.831337837423896,
    .area = 0.78,
  },
};
void test_raster(PolygonData const& testcase, Polygon const& poly, auto& test) {
    for (auto p : testcase.points) {
      test.check(poly.containsPoint(p));
    }

    if (!testcase.do_raster) {
      return;
    }

    std::string image = "(\n";
    for (size_t y_sample = testcase.samples_y; y_sample--; ) {
      std::string line;
      for (size_t x_sample = 0; x_sample < testcase.samples_x; ++x_sample) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
        double x = (testcase.x0 * (testcase.samples_x - x_sample - 1) + testcase.x1 * x_sample) / (testcase.samples_x - 1);
        double y = (testcase.y0 * (testcase.samples_y - y_sample - 1) + testcase.y1 * y_sample) / (testcase.samples_y - 1);
#pragma GCC diagnostic pop
        line += poly.containsPoint({x, y}) ? '#' : ' ';
      }
      while (line.size() && line.back() == ' ') {
        line.pop_back();
      }
      image += line + '\n';
    }
    image += ")";

    test.equals(testcase.image, image);
    if (testcase.image != image) {
      std::cout << "Expected image: " << testcase.image << "\nProduced image: " << image << std::endl;
    }
}
struct LinearTransform {
  double a0 = 1, b0 = 0, c0 = 0;
  double a1 = 0, b1 = 1, c1 = 0;

  bool is_identity = true;
  bool is_movement = true;

  Point applyTo(Point const& p) const {
    return {p.x * a0 + p.y * b0 + c0, p.x * a1 + p.y * b1 + c1};
  }

  LinearTransform rotate(double angle) const {
    double sin = std::sin(angle), cos = std::cos(angle);

    return {
      .a0 = a0 * cos - a1 * sin,
      .b0 = b0 * cos - b1 * sin,
      .c0 = c0 * cos - c1 * sin,

      .a1 = a0 * sin + a1 * cos,
      .b1 = b0 * sin + b1 * cos,
      .c1 = c0 * sin + c1 * cos,

      .is_identity = false,
      .is_movement = is_movement,
    };
  }

  LinearTransform move(Point p) const {
    return {
      .a0 = a0, .b0 = b0, .c0 = c0 + p.x,
      .a1 = a1, .b1 = b1, .c1 = c1 + p.y,

      .is_identity = false,
      .is_movement = is_movement,
    };
  }

  LinearTransform scale(double k) const {
    return {
      .a0 = k * a0, .b0 = k * b0, .c0 = k * c0,
      .a1 = k * a1, .b1 = k * b1, .c1 = k * c1,

      .is_identity = false,
      .is_movement = false,
    };
  }
};

constexpr LinearTransform IDENTITY{};
constexpr LinearTransform MIRROR{-1, 0, 0, 0, 1, 0, false};

const std::vector<LinearTransform> transforms_collection = {
  IDENTITY,
  IDENTITY.rotate(M_PI / 2),
  IDENTITY.move({1, 2}),
  IDENTITY.move({0, 1}),
  IDENTITY.scale(2),

  IDENTITY.move({-1, -1}).scale(2).move({1, 1}),
  IDENTITY.move({-1, 2}).rotate(M_PI / 3).move({1, -2}),

  MIRROR,
  MIRROR.rotate(M_PI / 2),
  MIRROR.move({1, 2}),
  MIRROR.move({0, 1}),
  MIRROR.scale(2),

  MIRROR.move({-1, -1}).scale(2).move({1, 1}),
  MIRROR.move({-1, 2}).rotate(M_PI / 3).move({1, -2}),

};

//NOLINTNEXTLINE
static TestGroup dan_tests[] = {

// ===== Point =====
{"Point",
make_test<PrettyTest>("ConstructionAndComparison", [](auto& test) {
  test.check((Point{1, 2} == Point{1, 2}));
  test.check(!(Point{1, 2} == Point{1.5, 2}));
  test.check(!(Point{1, 2} == Point{1, 2.5}));
  test.check(!(Point{1, 2} == Point{1.5, 2.5}));
  test.check((Point{0, 0} == Point{0, 0}));

  test.check(!(Point{1, 2} != Point{1, 2}));
  test.check((Point{1, 2} != Point{1.5, 2}));
  test.check((Point{1, 2} != Point{1, 2.5}));
  test.check((Point{1, 2} != Point{1.5, 2.5}));
  test.check(!(Point{0, 0} != Point{0, 0}));
})},

{"Line",
make_test<PrettyTest>("ConstructionAndComparison", [](auto& test) {
  std::vector<Line> lines = {
    /*0*/ Line{2, 5},
    /*1*/ Line{Point{0, 5}, 2},
    /*2*/ Line{Point{-4, -3}, 2},
    /*3*/ Line{Point{-4, -3}, Point{1.5, 8}},
    /*4*/ Line{Point{1.5, 8}, Point{-4, -3}},

    /*5*/ Line{2, 7},
    /*6*/ Line{Point{-0.5, 6}, 2},
    /*7*/ Line{Point{0, 7}, Point{-3.25, 0.5}},
    /*8*/ Line{Point{-3.25, 0.5}, 2},

    /*9*/ Line{0, 0},
    /*10*/ Line{Point{-0.3, 0}, Point{0.4, 0}},
    /*11*/ Line{Point{-224, 0}, 0},

    /*12*/ Line{0, 1},
    /*13*/ Line{Point{-0.3, 1}, Point{0.4, 1}},
    /*14*/ Line{Point{-224, 1}, 0},

    /*15*/ Line{Point{0, 0}, Point{0, 1}},
    /*16*/ Line{Point{0, 224}, Point{0, 37700}},

    /*17*/ Line{Point{1, 0}, Point{1, 1}},
    /*18*/ Line{Point{1, 420}, Point{1, 3665}},
  };

  std::vector<int> equiv_class = {
    0, 0, 0, 0, 0,
    1, 1, 1, 1,
    2, 2, 2,
    3, 3, 3,
    4, 4,
    5, 5,
  };

  for (size_t i = 0; i < lines.size(); ++i) {
    for (size_t j = 0; j < lines.size(); ++j) {
      test.equals(lines[i] == lines[j], equiv_class[i] == equiv_class[j]);
      test.equals(lines[i] != lines[j], equiv_class[i] != equiv_class[j]);
    }
  }
})},

{"Line",

make_test<PrettyTest>("Construction", [](auto& test){
  for (auto const& testcase : polygon_collection) {
    Polygon poly{testcase.points};
    test.equals(poly.verticesCount(), testcase.points.size());
    test.equals(poly.getVertices(), testcase.getPointsVector());
  }

  auto poly = Polygon{
      Point{0.8, 0.8},
      Point{2.6, 2.8},
      Point{0.6, 2.2},
  };
  test.equals(poly.verticesCount(), size_t(3));
  test.equals(poly.getVertices(), manual_test_points);
}),

make_test<PrettyTest>("GetVerticesOfRValue", [](auto& test){
  for (auto const& testcase : polygon_collection) {
    test.equals(Polygon{testcase.points}.getVertices(), testcase.getPointsVector());
  }

  test.equals((Polygon{
      Point{0.8, 0.8},
      Point{2.6, 2.8},
      Point{0.6, 2.2},
    }).getVertices(), manual_test_points);
}),

make_test<PrettyTest>("PerimeterAndArea", [](auto& test){
  for (auto const& testcase : polygon_collection) {
    test.check(std::abs(Polygon{testcase.points}.perimeter() - testcase.perimeter) < 1e-5);
    test.check(std::abs(Polygon{testcase.points}.area() - testcase.area) < 1e-5);
  }
}),

make_test<PrettyTest>("IsConvex", [](auto& test){
  for (auto const& testcase : polygon_collection) {
    test.equals(Polygon{testcase.points}.isConvex(), testcase.is_convex);
    test.equals(Polygon(testcase.getReversedPoints()).isConvex(), testcase.is_convex);
  }
}),
make_test<PrettyTest>("ContainsPoint", [](auto& test){
  for (auto const& testcase : polygon_collection) {
    Polygon poly{testcase.points};

    test_raster(testcase, poly, test);
    test_raster(testcase, Polygon(testcase.getReversedPoints()), test);
  }
}),
make_test<PrettyTest>("TrivialComparisons", [](auto& test){
  for (auto const& poly1 : polygon_collection) {
    for (auto const& poly2 : polygon_collection) {
      bool are_equal = &poly1 - polygon_collection == &poly2 - polygon_collection;

      test.equals(Polygon(poly1.points) == Polygon(poly2.points), are_equal);
      test.equals(Polygon(poly1.points).isCongruentTo(Polygon(poly2.points)), are_equal);
      test.equals(Polygon(poly1.points).isSimilarTo(Polygon(poly2.points)), are_equal);
    }
  }
}),
make_test<PrettyTest>("Comparisons", [](auto& test){
  for (auto const& poly : polygon_collection) {
    for (auto const& poly1 : polygon_collection) {
      bool are_equal = &poly - polygon_collection == &poly1 - polygon_collection;

      for (auto const& transform : transforms_collection) {
        std::vector<Point> poly2 = poly1.getPointsVector();
        for (auto& p : poly2) {
          p = transform.applyTo(p);
        }

        for (size_t is_reversed = 0; is_reversed < 2; ++is_reversed) {
          for (size_t shift = 0; shift < poly2.size(); ++shift) {
            Polygon poly3(poly.points);
            test.equals(poly3 == Polygon(poly2), transform.is_identity && are_equal);
            test.equals(poly3.isCongruentTo(Polygon(poly2)), transform.is_movement && are_equal);
            test.equals(poly3.isSimilarTo(Polygon(poly2)), are_equal);
            std::rotate(poly2.begin(), poly2.begin() + 1, poly2.end());
          }
          std::reverse(poly2.begin(), poly2.end());
        }
      }
    }
  }
}),

make_test<PrettyTest>("Transforms", [](auto& test){
  for (auto const& poly : polygon_collection) {
    auto check = [&](Polygon const& got, LinearTransform expected_transform, std::string what) {
      std::vector<Point> short_cut = poly.getPointsVector();
      for (auto& p : short_cut) {
        p = expected_transform.applyTo(p);
      }

      if (!test.equals(got, Polygon{short_cut})) {
        std::cout << what << '\n';
      }
    };

    Polygon p(poly.points), q = p;
    q.rotate({0, 0}, 45);
    check(q, IDENTITY.rotate(M_PI / 4), "rotate 45deg");

    q = p;
    q.rotate({1, 2}, -45);
    check(q, IDENTITY.move({-1, -2}).rotate(-M_PI / 4).move({1, 2}), "rotate 45deg around (1, 2)");

    q = p;
    q.reflect(Point{0, 0});
    check(q, {-1, 0, 0, 0, -1, 0}, "reflect over (0, 0)");

    q = p;
    q.reflect(Point{1, 2});
    check(q, {-1, 0, 2, 0, -1, 4}, "reflect over (1, 2)");

    q = p;
    q.scale({0, 0}, 2);
    check(q, {2, 0, 0, 0, 2, 0}, "scale by 2");

    q = p;
    q.scale({1, 2}, 2);
    check(q, {2, 0, -1, 0, 2, -2}, "scale by 2 using (1, 2) as the origin");

    q = p;
    q.reflect(Line{0, 0});
    check(q, {1, 0, 0, 0, -1, 0}, "reflect through y=0");

    q = p;
    q.reflect(Line{Point{0, 0}, Point{0, 2}});
    check(q, {-1, 0, 0, 0, 1, 0}, "reflect through x=0");

    q = p;
    q.reflect(Line{Point{1, 0}, Point{0, 1}});
    check(q, {0, -1, 1, -1, 0, 1}, "reflect through x+y=1 (#1)");

    q = p;
    q.reflect(Line{Point{2, -1}, Point{0, 1}});
    check(q, {0, -1, 1, -1, 0, 1}, "reflect through x+y=1 (#2)");
  }
})
},

{"Ellipse",
make_test<PrettyTest>("ConstructionAndGetters", [](auto& test){
  Ellipse a{{-1.5, 0}, {1.5, 0}, 5};
  Ellipse b{{1, 1}, {4, 5}, 6};

  test.equals(a.focuses(), (std::pair<Point, Point>{{-1.5, 0}, {1.5, 0}}));
  test.equals(a.directrices(), (std::pair<Line, Line>{
    {{-4.166666666666667, 0}, {-4.166666666666667, 1}},
    {{4.166666666666667, 0}, {4.166666666666667, 1}}
  }));
  test.equals(a.center(), (Point{0, 0}));
  test.float_equals(a.eccentricity(), 0.6, double_error);
  test.float_equals(a.area(), 15.707963267948966, double_error);
  {
      const double real_perimeter = 14.180833944487244;
      const double perimeter_error = real_perimeter * 0.05 + double_error;
      test.float_equals(a.perimeter(), real_perimeter, perimeter_error);
  }
  test.equals(b.focuses(), (std::pair<Point, Point>{{1, 1}, {4, 5}}));
  test.equals(b.directrices(), (std::pair{
    Line{-0.75, 0.375},
    Line{-0.75, 9.375},
  }));
  test.equals(b.center(), (Point{2.5, 3}));
  test.float_equals(b.eccentricity(), 5. / 6, double_error);
  test.float_equals(b.area(), 15.629226114141467, double_error);
  {
      const double real_perimeter = 14.180833944487244;
      const double perimeter_error = real_perimeter * 0.05 + double_error;
      test.float_equals(b.perimeter(), 14.939631086609392, perimeter_error);
  }
})
},

{"Rectangle",
make_test<PrettyTest>("Construction", [](auto& test) {
  test.equals((Rectangle{{1, 1}, {3, 3}, 1}), (Polygon{{1, 1}, {3, 1}, {3, 3}, {1, 3}}));
  test.equals((Rectangle{{1, 2}, {5, 5}, 0.75}), (Polygon{{1, 2}, {5, 2}, {5, 5}, {1, 5}}));
  test.equals((Rectangle{{1, 2}, {5, 5}, 4 / 3.}), (Polygon{{1, 2}, {5, 2}, {5, 5}, {1, 5}}));
})

},

{"Triangle",
make_test<PrettyTest>("Construction and getters", [](auto& test) {
  Triangle a{{3, 0}, {0, 3}, {6, 1}};

  test.equals(a.circumscribedCircle(), (Circle{{3.5, 3.5}, 3.5355339059327373}));
  test.equals(a.inscribedCircle(), (Circle{{3.23606797749979, 1}, 0.8740320488976421}));
  test.equals(a.centroid(), (Point{3, 1.3333333333333333}));
  test.equals(a.orthocenter(), (Point{2, -3}));
  test.equals(a.EulerLine(), (Line{{2, -3}, {3, 1.3333333333333333}}));
  test.equals(a.ninePointsCircle(), (Circle{{2.75, 0.25}, 1.7677669529663687}));

})
}

};

