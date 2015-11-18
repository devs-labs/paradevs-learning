/**
 * @file tests/pdevs/graph_manager.cpp
 * @author The PARADEVS Development Team
 * See the AUTHORS or Authors.txt file
 */

/*
 * PARADEVS - the multimodeling and simulation environment
 * This file is a part of the PARADEVS environment
 *
 * Copyright (C) 2013-2015 ULCO http://www.univ-litoral.fr
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TESTS_PDEVS_GRAPH_MANAGER_HPP
#define TESTS_PDEVS_GRAPH_MANAGER_HPP 1

#include <models.hpp>

#include <paradevs/kernel/pdevs/Coordinator.hpp>
#include <paradevs/kernel/pdevs/GraphManager.hpp>
#include <paradevs/kernel/pdevs/Simulator.hpp>

namespace paradevs { namespace tests { namespace pdevs {

template < typename M >
class OnlyOneGraphManager :
        public paradevs::pdevs::GraphManager < common::DoubleTime >
{
public:
    OnlyOneGraphManager(common::Coordinator < common::DoubleTime >* coordinator,
                        const paradevs::common::NoParameters& parameters) :
        paradevs::pdevs::GraphManager < common::DoubleTime >(coordinator,
                                                          parameters),
        model("m", common::NoParameters())
    {
        add_child(&model);
    }

    virtual ~OnlyOneGraphManager()
    { }

private:
    paradevs::pdevs::Simulator < common::DoubleTime, M > model;
};

class TwoModelsGraphManager :
        public paradevs::pdevs::GraphManager < common::DoubleTime >
{
public:
    TwoModelsGraphManager(common::Coordinator < common::DoubleTime >*
                          coordinator,
                          const paradevs::common::NoParameters& parameters) :
        paradevs::pdevs::GraphManager < common::DoubleTime >(coordinator,
                                                             parameters),
        a("a", common::NoParameters())
        // ,
        // b("b", common::NoParameters()),
        // c("c", common::NoParameters()),
        // d("d", common::NoParameters()),
        // e("e", common::NoParameters()),
        // f("f", common::NoParameters()),
        // g("g", common::NoParameters()),
        // h("h", common::NoParameters())
    {
        add_child(&a);
        // add_child(&b);
        // add_child(&c);
        // add_child(&d);
        // add_child(&e);
        // add_child(&f);
        // add_child(&g);
        // add_child(&h);

        // a.add_out_port("out");
        // b.add_in_port("in");
        // b.add_out_port("out");
        // c.add_in_port("in");
        // c.add_out_port("out");
        // d.add_in_port("in");
        // d.add_out_port("out");
        // e.add_in_port("in");
        // e.add_out_port("out");
        // f.add_in_port("in");
        // f.add_out_port("out");
        // g.add_in_port("in");
        // g.add_out_port("out");
        // h.add_in_port("in");
        // h.add_out_port("out");

        // add_link(&a, "out", &d, "in");
        // add_link(&b, "out", &d, "in");
        // add_link(&b, "out", &e, "in");
        // add_link(&c, "out", &e, "in");
        // add_link(&d, "out", &f, "in");
        // add_link(&e, "out", &f, "in");
        // add_link(&e, "out", &g, "in");
        // add_link(&f, "out", &h, "in");
        // add_link(&g, "out", &h, "in");
    }

    virtual ~TwoModelsGraphManager()
    { }

private:
    paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > a;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > b;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > c;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > d;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > e;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > f;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > g;
    // paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > h;
};

class GeneratorGraphManager :
        public paradevs::pdevs::GraphManager < common::DoubleTime >
{
public:
    GeneratorGraphManager(common::Coordinator < common::DoubleTime >*
                          coordinator,
                          const paradevs::common::NoParameters& parameters) :
        paradevs::pdevs::GraphManager < common::DoubleTime >(coordinator,
                                                             parameters),
        a("a", common::NoParameters()), b("b", common::NoParameters())
    {
        add_child(&a);
        add_child(&b);

        a.add_out_port("out");
        b.add_in_port("in");
        b.add_out_port("out");

        add_link(&a, "out", &b, "in");
    }

    virtual ~GeneratorGraphManager()
    { }

private:
    paradevs::pdevs::Simulator < common::DoubleTime, Generator > a;
    paradevs::pdevs::Simulator < common::DoubleTime, ThreeStateModel > b;
};

} } } // namespace paradevs tests pdevs

#endif
