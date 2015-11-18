/**
 * @file tests/pdevs/models.hpp
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

#ifndef TESTS_PDEVS_MODELS_HPP
#define TESTS_PDEVS_MODELS_HPP 1

#include <paradevs/common/time/DoubleTime.hpp>

#include <paradevs/kernel/pdevs/Dynamics.hpp>

namespace paradevs { namespace tests { namespace pdevs {

class TwoStateModel :
        public paradevs::pdevs::Dynamics < common::DoubleTime >
{
public:
    TwoStateModel(const std::string& name,
                  const common::NoParameters& parameters) :
        paradevs::pdevs::Dynamics < common::DoubleTime >(name, parameters)
    { }
    virtual ~TwoStateModel()
    { }

    void dint(typename common::DoubleTime::type t)
    {
        if (_phase == S1) {
            _phase = S2;
        } else if (_phase == S2) {
            _phase = S1;
        }
        _last_time = t;
    }

    typename common::DoubleTime::type start(
        typename common::DoubleTime::type t)
    {
        _phase = S1;
        _last_time = t;
        return ta(t);
    }

    typename common::DoubleTime::type ta(
        typename common::DoubleTime::type /* t */) const
    {
        if (_phase == S1) {
            return 5;
        } else {
            return 6;
        }
    }

    common::Bag < common::DoubleTime > lambda(
        typename common::DoubleTime::type t) const
    {

        std::cout << (t - _last_time) << std::endl;

        return common::Bag < common::DoubleTime >();
    }

private:
    enum Phase { S1, S2 };

    Phase _phase;
    typename common::DoubleTime::type _last_time;
};

class ThreeStateModel :
        public paradevs::pdevs::Dynamics < common::DoubleTime >
{
public:
    ThreeStateModel(const std::string& name,
                    const common::NoParameters& parameters) :
        paradevs::pdevs::Dynamics < common::DoubleTime >(name, parameters)
    { }
    virtual ~ThreeStateModel()
    { }

    void compute()
    {
        for (unsigned int i = 0; i < heights.size(); ++i) {
            if (heights[i] != -1 and heights[i] < 10) {
                heights[i] += speeds[i] * scales[i];
            }
        }
    }

    void display() const
    {
        for (std::vector < double >::const_iterator it = heights.begin();
             it != heights.end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout << std::endl;
    }

    void display_full() const
    {
        unsigned int i = 1;

        for (std::vector < double >::const_iterator it = heights.begin();
             it != heights.end(); ++it, ++i) {
            if (*it > 10) {
                std::cout << "S" << i;
            }
        }
        std::cout << std::endl;
    }

    bool full() const
    {
        unsigned int n = 0;

        for (std::vector < double >::const_iterator it = heights.begin();
             it != heights.end(); ++it) {
            if (*it > 10) {
                ++n;
            }
        }
        return n > 0;
    }

    bool full_N() const
    {
        unsigned int n = 0;

        for (std::vector < double >::const_iterator it = heights.begin();
             it != heights.end(); ++it) {
            if (*it == -1) {
                ++n;
            }
        }
        return n >= 2;
    }

    void mark_full(typename common::DoubleTime::type t)
    {
        for (std::vector < double >::iterator it = heights.begin();
             it != heights.end(); ++it) {
            if (*it > 10) {
                *it = -1;
                _last_time = t;
            }
        }
    }

    void raz()
    {
        for (std::vector < double >::iterator it = heights.begin();
             it != heights.end(); ++it) {
            if (*it == -1) {
                *it = 0;
            }
        }
    }

    void dconf(typename common::DoubleTime::type t,
               typename common::DoubleTime::type e,
               const common::Bag < common::DoubleTime >& msgs)
    {
        dext(t, e, msgs);
    }

    void dext(typename common::DoubleTime::type t,
              typename common::DoubleTime::type /* e */,
              const common::Bag < common::DoubleTime >& msgs)
    {
        for (common::Bag < common::DoubleTime >::const_iterator
                 it = msgs.begin(); it != msgs.end(); ++it) {
            ++n;
        }

        // std::cout << get_name() << "|I|" << (t - _last_time) << " ";

        if (sigma == 1) {
            if (n > 3) {
                ++index;
                if (index == scales.size()) {
                    index = 0;
                }
                sigma = std::numeric_limits < double >::max();
                if (scales[index] == 1) {
                    scales[index] = 2;
                } else {
                    scales[index] = 1;
                }

                std::cout << get_name() << "||INF ";

                n = 0;
            }
        } else {
            sigma = 1;
            n = 0;
        }
    }

    void dint(typename common::DoubleTime::type t)
    {
        mark_full(t);
        if (full_N()) {
            raz();
        }
        compute();
    }

    typename common::DoubleTime::type start(
        typename common::DoubleTime::type t)
    {
        heights = { 0, 0, 0, 0, 0 };
        speeds = { 0.21, 0.3, 0.7, 0.56, 0.14 };
        scales = { 1, 1, 1, 1, 1 };
        index = 0;
        n = 0;
        sigma = 1;
        _last_time = t;
        return 0;
    }

    typename common::DoubleTime::type ta(
        typename common::DoubleTime::type /* t */) const
    { return sigma; }

    common::Bag < common::DoubleTime > lambda(
        typename common::DoubleTime::type t) const
    {
        common::Bag < common::DoubleTime > msgs;

        if (full()) {

            std::cout << get_name() << "|O|" << (t - _last_time) << " ";
            // display_full();

            msgs.push_back(common::ExternalEvent < common::DoubleTime >(
                               "out", 0));
        }
        return msgs;
    }

private:
    std::vector < double > heights;
    std::vector < double > speeds;
    std::vector < double > scales;
    unsigned int index;
    unsigned int n;
    typename common::DoubleTime::type sigma;

    typename common::DoubleTime::type _last_time;
};

class Generator :
        public paradevs::pdevs::Dynamics < common::DoubleTime >
{
public:
    Generator(const std::string& name,
              const common::NoParameters& parameters) :
        paradevs::pdevs::Dynamics < common::DoubleTime >(name, parameters)
    { }
    virtual ~Generator()
    { }

    typename common::DoubleTime::type start(
        typename common::DoubleTime::type /* t */)
    { return rand() % 20 + 1; }

    typename common::DoubleTime::type ta(
        typename common::DoubleTime::type /* t */) const
    { return rand() % 20 + 1; }

    common::Bag < common::DoubleTime > lambda(
        typename common::DoubleTime::type /* t */) const
    {
        common::Bag < common::DoubleTime > msgs;

        msgs.push_back(common::ExternalEvent < common::DoubleTime >(
                           "out", 0));
        return msgs;
    }
};

} } } // namespace paradevs tests pdevs

#endif
