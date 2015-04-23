//============================================================================
// test_random_access.cpp
// Copyright (c) 2011 Pete Goodliffe. All rights reserved
//============================================================================

#include "random_access_skip_list.h"

#define CATCH_CONFIG_NO_STREAM_REDIRECTION 1
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "test_types.h"

using goodliffe::random_access_skip_list;

TEST_CASE( "random_access_skip_list/smoketest", "" )
{
    //REQUIRE(false);
}

TEST_CASE( "random_access_skip_list/is constructable", "" )
{
    random_access_skip_list<int>    sl_int;
    random_access_skip_list<float>  sl_float;
    random_access_skip_list<Struct> sl_struct;
    
    REQUIRE(true);
}

TEST_CASE( "random_access_skip_list/can call basic methods from skip_list", "" )
{
    const random_access_skip_list<int> list;
    REQUIRE(list.empty());
    REQUIRE(list.size() == 0);
    REQUIRE(list.get_allocator() == std::allocator<int>());
    REQUIRE(list.find(10) == list.end());
    REQUIRE(list.count(0) == 0);
    REQUIRE(!list.contains(20));
    REQUIRE(list.max_size() > 67890);
    REQUIRE(list.begin() == list.end());
    REQUIRE(list.rbegin() == list.rend());
}

//============================================================================
// indexing

TEST_CASE( "random_access_skip_list/add items", "" )
{
    random_access_skip_list<int> list;
    REQUIRE(list.size() == 0);

    list.insert(0);
    REQUIRE(list.size() == 1);
    REQUIRE(list[0] == 0);

    list.insert(1);
    REQUIRE(list.size() == 2);
    REQUIRE(list[0] == 0);
    REQUIRE(list[1] == 1);
    
    list.insert(2);
    REQUIRE(list.size() == 3);
    REQUIRE(list[0] == 0);
    REQUIRE(list[1] == 1);
    REQUIRE(list[2] == 2);
    
    list.insert(3);
    REQUIRE(list.size() == 4);
    REQUIRE(list[0] == 0);
    REQUIRE(list[1] == 1);
    REQUIRE(list[2] == 2);
    REQUIRE(list[3] == 3);
    
    list.insert(4);
    REQUIRE(list.size() == 5);
    REQUIRE(list[0] == 0);
    REQUIRE(list[1] == 1);
    REQUIRE(list[2] == 2);
    REQUIRE(list[3] == 3);
    REQUIRE(list[4] == 4);
    
    list.insert(5);
    REQUIRE(list.size() == 6);
    REQUIRE(list[0] == 0);
    REQUIRE(list[1] == 1);
    REQUIRE(list[2] == 2);
    REQUIRE(list[3] == 3);
    REQUIRE(list[4] == 4);
}

TEST_CASE( "random_access_skip_list/indexing", "" )
{
    for (unsigned n = 0; n < 10; ++n)
    {
        unsigned size = unsigned(rand()) % 300;
        std::vector<int> data;
        random_access_skip_list<int> list;
        FillWithRandomData(size, data);
        for (unsigned m = 0 ; m < size; ++m)
            list.insert(data[m]);
        
        SortVectorAndRemoveDuplicates(data);
        
        for (unsigned m = 0 ; m < data.size(); ++m)
        {
            REQUIRE(list[m] == data[m]);
        }

        // Check that clear puts indexing back to right state so subsequent
        // insertions are correct
        list.clear();
        REQUIRE(list.size() == 0);
        REQUIRE(list.empty());

        for (unsigned m = 0 ; m < data.size(); ++m)
            list.insert(data[data.size()-1-m]);
        
        for (unsigned m = 0 ; m < data.size(); ++m)
        {
            REQUIRE(list[m] == data[m]);
        }
    }
}

//============================================================================
#pragma mark iterators

TEST_CASE( "random_access_skip_list/iterators/can be default constructed", "" )
{
    random_access_skip_list<int>::iterator i;
    random_access_skip_list<int>::const_iterator ci;
}

TEST_CASE( "random_access_skip_list/iterators/convert from one to t'other", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);

    random_access_skip_list<int>::iterator i = list.begin();
    random_access_skip_list<int>::const_iterator ci = i;
}

TEST_CASE( "random_access_skip_list/iterator_at", "" )
{
    std::vector<int> data;
    FillWithRandomData(1000, data);
    
    random_access_skip_list<int> list(data.begin(), data.end());
    const random_access_skip_list<int> clist = list;
    
    SortVectorAndRemoveDuplicates(data);
    
    for (unsigned n = 0; n < list.size(); n += 7)
    {
        REQUIRE(*list.iterator_at(n) == data[n]);
        REQUIRE(*list.citerator_at(n) == data[n]);
        REQUIRE(*clist.iterator_at(n) == data[n]);
    }
}

TEST_CASE( "random_access_skip_list/iterators/+=", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    for (int n = 0; n < 9; ++n) list.insert(n);
    
    random_access_skip_list<int>::iterator i = list.begin();
    REQUIRE(*i == 0);
    i += 2;    REQUIRE(*i == 2);
    i += 2;    REQUIRE(*i == 4);
    i += 3;    REQUIRE(*i == 7);
    i += 1;    REQUIRE(*i == 8);
    
    random_access_skip_list<int>::const_iterator ci = clist.begin();
    REQUIRE(*ci == 0);
    ci += 2;    REQUIRE(*ci == 2);
    ci += 2;    REQUIRE(*ci == 4);
    ci += 3;    REQUIRE(*ci == 7);
    ci += 1;    REQUIRE(*ci == 8);
}

TEST_CASE( "random_access_skip_list/iterators/-=", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    for (int n = 0; n < 10; ++n) list.insert(n);

    REQUIRE(list.index_of(list.end()) == 10);
    REQUIRE(clist.index_of(list.end()) == 10);

    random_access_skip_list<int>::iterator i = list.end();
    i -= 2;    REQUIRE(*i == 8);
    i -= 2;    REQUIRE(*i == 6);
    i -= 3;    REQUIRE(*i == 3);
    i -= 1;    REQUIRE(*i == 2);

    random_access_skip_list<int>::const_iterator ci = clist.end();
    ci -= 2;    REQUIRE(*ci == 8);
    ci -= 2;    REQUIRE(*ci == 6);
    ci -= 3;    REQUIRE(*ci == 3);
    ci -= 1;    REQUIRE(*ci == 2);
}

TEST_CASE( "random_access_skip_list/iterators/+/lhs is iterator", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    std::vector<int> vec;
    for (int n = 0; n < 10; ++n) list.insert(n);
    for (int n = 0; n < 10; ++n) vec.push_back(n);
    
    {
        // A comparison of what can be done with vector
        std::vector<int>::iterator i = vec.begin() + 2;
        REQUIRE(*i == 2);
    }
    {
        random_access_skip_list<int>::iterator i = list.begin() + 2;
        REQUIRE(*i == 2);
        random_access_skip_list<int>::const_iterator ci = clist.begin() + 2;
        REQUIRE(*ci == 2);
    }
}

TEST_CASE( "random_access_skip_list/iterators/+/rhs is iterator", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    std::vector<int> vec;
    for (int n = 0; n < 10; ++n) list.insert(n);
    for (int n = 0; n < 10; ++n) vec.push_back(n);
    
    {
        // A comparison of what can be done with vector
        std::vector<int>::iterator i = 2 + vec.begin();
        REQUIRE(*i == 2);
    }
    {
        random_access_skip_list<int>::iterator i = 2 + list.begin();
        REQUIRE(*i == 2);
        random_access_skip_list<int>::const_iterator ci = 2 + clist.begin();
        REQUIRE(*ci == 2);
    }
}

TEST_CASE( "random_access_skip_list/iterators/-/lhs is iterator", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    std::vector<int> vec;
    for (int n = 0; n < 10; ++n) list.insert(n);
    for (int n = 0; n < 10; ++n) vec.push_back(n);
    
    {
        // A comparison of what can be done with vector
        std::vector<int>::iterator i = vec.begin() + 5;
        i = i - 2;
        REQUIRE(*i == 3);
    }
    {
        random_access_skip_list<int>::iterator i = list.begin() + 5;
        i = i - 2;
        REQUIRE(*i == 3);
    }
    {
        random_access_skip_list<int>::const_iterator i = clist.begin() + 5;
        i = i - 2;
        REQUIRE(*i == 3);
    }
}

TEST_CASE( "random_access_skip_list/iterators/operator[]", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    std::vector<int> vec;
    for (int n = 0; n < 10; ++n) list.insert(n);
    for (int n = 0; n < 10; ++n) vec.push_back(n);
    
    {
        // A comparison of what can be done with vector
        std::vector<int>::iterator i = vec.begin()+2;
        REQUIRE(i[3] == 5);
    }
    {
        random_access_skip_list<int>::iterator i = list.begin()+2;
        REQUIRE(i[3] == 5);
    }
    {
        random_access_skip_list<int>::const_iterator i = clist.begin()+2;
        REQUIRE(i[3] == 5);
    }
}

TEST_CASE( "random_access_skip_list/iterators/operator<", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    std::vector<int> vec;
    for (int n = 0; n < 10; ++n) list.insert(n);
    for (int n = 0; n < 10; ++n) vec.push_back(n);
    
    {
        // A comparison of what can be done with vector
        REQUIRE(vec.begin() < vec.end());
        REQUIRE_FALSE(vec.end() < vec.begin());
        
        REQUIRE_FALSE(vec.begin() < vec.begin());
        REQUIRE(vec.begin() < vec.begin()+1);
    }
    {
        REQUIRE(list.begin() < list.end());
        REQUIRE_FALSE(list.end() < list.begin());
        
        REQUIRE_FALSE(list.begin() < list.begin());
        REQUIRE(list.begin() < list.begin()+1);
    }
    {
        REQUIRE(clist.begin() < clist.end());
        REQUIRE_FALSE(clist.end() < clist.begin());
        
        REQUIRE_FALSE(clist.begin() < clist.begin());
        REQUIRE(clist.begin() < clist.begin()+1);
    }
}

TEST_CASE( "random_access_skip_list/iterators/subtraction", "" )
{
    random_access_skip_list<int> list;
    const random_access_skip_list<int> &clist = list;
    for (int n = 0; n < 10; ++n) list.insert(n);
    
    REQUIRE((list.end()-list.begin()) == 10);
    REQUIRE(((list.begin()+5)-(list.begin()+2)) == 3);

    REQUIRE((clist.end()-clist.begin()) == 10);
    REQUIRE(((clist.begin()+5)-(clist.begin()+2)) == 3);
}

//============================================================================
#pragma mark erase_at

TEST_CASE( "random_access_skip_list/erase_at/maintains indexes/1", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);
    
    list.erase_at(1);
    REQUIRE(list[0] == 0);    REQUIRE(list[1] == 2);    REQUIRE(list[2] == 3);
    REQUIRE(list[3] == 4);    REQUIRE(list[4] == 5);    REQUIRE(list[7] == 8);
}

TEST_CASE( "random_access_skip_list/erase_at/maintains indexes/2", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);
    
    list.erase_at(2);
    REQUIRE(list[0] == 0);    REQUIRE(list[1] == 1);    REQUIRE(list[2] == 3);
    REQUIRE(list[3] == 4);    REQUIRE(list[4] == 5);    REQUIRE(list[7] == 8);
}

TEST_CASE( "random_access_skip_list/erase_at/maintains indexes/3", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);

    list.erase_at(3);
    REQUIRE(list[0] == 0);    REQUIRE(list[1] == 1);    REQUIRE(list[2] == 2);
    REQUIRE(list[3] == 4);    REQUIRE(list[4] == 5);    REQUIRE(list[7] == 8);
}

TEST_CASE( "random_access_skip_list/erase_at/maintains indexes/4", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);
    
    list.erase_at(4);
    REQUIRE(list[0] == 0);    REQUIRE(list[1] == 1);    REQUIRE(list[2] == 2);
    REQUIRE(list[3] == 3);    REQUIRE(list[4] == 5);    REQUIRE(list[7] == 8);
}

TEST_CASE( "random_access_skip_list/erase_at/maintains indexes", "" )
{
    std::vector<int> data;
    FillWithRandomData(1000, data);
    
    random_access_skip_list<int> list(data.begin(), data.end());

    SortVectorAndRemoveDuplicates(data);

    for (unsigned n = 0 ; n < 20; ++n)
    {
        unsigned index = unsigned(rand()) % unsigned(list.size());
        //fprintf(stderr, "  remove %u/%u\n", index, (unsigned)list.size());
        list.erase_at(index);
        data.erase(data.begin()+index);
    }
    
    REQUIRE(CheckForwardIteration(list));
    REQUIRE(CheckBackwardIteration(list));
    
    REQUIRE(list.size() == data.size());
    /*
    for (unsigned n = 0; n < list.size(); ++n)
    {
        fprintf(stderr, "  %u/%u: %d==%d %s\n", n, (unsigned)list.size(), list[n], data[n], list[n] == data[n]?"":"NO");
    }
     */
    REQUIRE(CheckEquality(list, data));
}

//============================================================================

TEST_CASE( "random_access_skip_list/erase range/maintains indexes", "" )
{
    for (unsigned megaloop = 1; megaloop < 10; ++megaloop)
    {
        std::vector<int> data;
        FillWithOrderedData(10 * megaloop, data);
        random_access_skip_list<int> list(data.begin(), data.end());

        for (unsigned m = 0; m < 5; ++m)
        {
            //list.dump(std::cerr);
            unsigned start  = unsigned(rand()) % unsigned(data.size()/2);
            unsigned length = unsigned(rand()) % unsigned(data.size()/2);
            //fprintf(stderr, "removing: (%u,%u)...\n", start, length);
            
            random_access_skip_list<int>::iterator i1 = list.iterator_at(start);
            random_access_skip_list<int>::iterator i2 = list.iterator_at(start+length);
            list.erase(i1, i2);
            
            std::vector<int>::iterator v1 = data.begin()+start;
            std::vector<int>::iterator v2 = data.begin()+start+length;
            data.erase(v1, v2);

            REQUIRE(CheckEquality(list, data));

            /*
            list.dump(std::cerr);
            for (unsigned n = 0; n < data.size(); ++n)
            {
                fprintf(stderr, "  %03u: d=% 6d/l=% 6d %s\n", n, data[n], list[n],
                        data[n] != list[n] ? "(error)":"");
            }
             */

            REQUIRE(CheckEqualityViaIndexing(list, data));
        }
    }
}

TEST_CASE( "random_access_skip_list/non members", "" )
{
    random_access_skip_list<int> list1, list2;
    list1.insert(1);
    list2.insert(2);
    
    REQUIRE(list1 != list2);
    REQUIRE_FALSE(list1 == list2);
    
    REQUIRE(list1 < list2);
    REQUIRE(list1 <= list2);
    REQUIRE_FALSE(list1 > list2);
    REQUIRE_FALSE(list1 >= list2);
}

//============================================================================
#pragma mark index_of

TEST_CASE( "random_access_skip_list/index_of", "" )
{
    random_access_skip_list<int> list;
    for (int n = 0; n < 9; ++n) list.insert(n);

    random_access_skip_list<int>::iterator i = list.begin();
    
    REQUIRE(list.index_of(i++) == 0);
    REQUIRE(list.index_of(i++) == 1);
    REQUIRE(list.index_of(i++) == 2);
    REQUIRE(list.index_of(i++) == 3);
    REQUIRE(list.index_of(i++) == 4);
    REQUIRE(list.index_of(i++) == 5);
    REQUIRE(list.index_of(i++) == 6);
    REQUIRE(list.index_of(i++) == 7);
    REQUIRE(list.index_of(i++) == 8);
    
    REQUIRE(list.index_of(list.end()) == 9);
}

//============================================================================
// TODO: allocation test

