//
//  TestMain.cpp
//
//  (C) Copyright Szymon Wojciechowski 2012-2013.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <iostream>
#include <utility>
#include <gtest/gtest.h>
//#include <gmock/gmock.h>
#define private public
#include <cosi/order_statistics.hpp>

//#undef private

#include <set>
#include <vector>

using cosi::util::order_statistics_tree;

typedef order_statistics_tree<int>::tree_node node_int;

template<typename ValueType, typename Comparator = std::less<ValueType> >
struct IteratorCmp
{
	typedef typename order_statistics_tree<ValueType, Comparator>::iterator iterator;
	bool operator()(const iterator& lhs, const iterator& rhs)
	{
		return Comparator()(*lhs, *rhs);
	}
};

struct LengthyStringComparator
{

	 typedef std::string first_argument_type;
	 typedef std::string second_argument_type;
	 typedef bool result_type;
	 
	bool operator()(const std::string& lhs, const std::string& rhs) const
	{
		return lhs.size() < rhs.size() || ( lhs.size() == rhs.size() && lhs < rhs );
	}
};

struct NodeCmp
{
	bool operator()(const node_int * a, const node_int * b)
	{
		return a->m_value < b->m_value;
	}
};

TEST(Insertion, RootInsertion)
{
	order_statistics_tree<int> t;
	node_int* root = new node_int(15);
	t.appendNode(root);
	ASSERT_TRUE(t.m_root == root);
	ASSERT_TRUE(root->m_parent == &t.m_end);
	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(1), t.size());
}

TEST(Insertion, LeftInsertion)
{
	order_statistics_tree<int> t;
	node_int* root = new node_int(15);
	node_int* node = new node_int(1);
	t.appendNode(root);
	t.appendNode(node);
	ASSERT_TRUE(t.m_root->m_left == node);
	ASSERT_TRUE(node->m_parent == root);
}

TEST(Insertion, RightInsertion)
{
	order_statistics_tree<int> t;
	node_int* root = new node_int(1);
	node_int* node = new node_int(15);
	t.appendNode(root);
	t.appendNode(node);
	ASSERT_TRUE(t.m_root->m_right == node);
	ASSERT_TRUE(node->m_parent == root);
}

TEST(Insertion, Insert)
{
	order_statistics_tree<int> t;
	int root = 15;
	order_statistics_tree<int>::iterator it = t.insert(root).first;
	ASSERT_EQ(root, *it);
	ASSERT_EQ(static_cast<size_t>(0), it.position());
}

TEST(Rotation, Left)
{
	order_statistics_tree<int> t;
	int values[] = {5, 15, 3, 18, 10};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
		t.appendNode(nodes[i]);
	}
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_left == nodes[4]);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[1]);
	
	t.leftRotation(nodes[0]);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[0]);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(5), t.size());
}

TEST(Rotation, Right)
{
	order_statistics_tree<int> t;
	int values[] = {15, 5, 3, 18, 10};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
		t.appendNode(nodes[i]);
	}

	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[1]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[1]);
	
	t.rightRotation(nodes[0]);

	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[1]->m_right == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == nodes[4]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[0]);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(5), t.size());
}

TEST(RBInsertion, In123)
{
	order_statistics_tree<int> t;
	int values[] = {1, 2, 3};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]); //case 3
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_isRed == true);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
}

TEST(RBInsertion, In132)
{
	order_statistics_tree<int> t;
	int values[] = {1, 3, 2};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]); //case 2,3
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
}

TEST(RBInsertion, In321)
{
	order_statistics_tree<int> t;
	int values[] = {3, 2, 1};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]); //case 3
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == true);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
}

TEST(RBInsertion, In312)
{
	order_statistics_tree<int> t;
	int values[] = {3, 1, 2};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]); //case 2,3
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_right == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
}

TEST(RBInsertion, In213)
{
	order_statistics_tree<int> t;
	int values[] = {2, 1, 3};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
}

TEST(RBInsertion, In231)
{
	order_statistics_tree<int> t;
	int values[] = {2, 3, 1};
	node_int* nodes[3];
	for(int i = 0 ; i < 3 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[0]->m_right == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);

	t.insertAndRebalance(nodes[2]);
	ASSERT_TRUE(t.m_root == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_right == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
}

TEST(RBInsertion, In1203)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3};
	node_int* nodes[4];
	for(int i = 0 ; i < 4 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == 0);

	t.insertAndRebalance(nodes[2]);
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_left == 0);

	t.insertAndRebalance(nodes[0]);
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_isRed == true);

	t.insertAndRebalance(nodes[3]); // case 1
	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[3]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[3]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
}

TEST(RBInsertion, In2310)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3};
	node_int* nodes[4];
	for(int i = 0 ; i < 4 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.insertAndRebalance(nodes[2]);
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_right == 0);

	t.insertAndRebalance(nodes[3]);
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[3]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_left == 0);

	t.insertAndRebalance(nodes[1]);
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);

	t.insertAndRebalance(nodes[0]); //case 1
	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[3]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[3]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[0]->m_isRed == true);
}

TEST(Neighbours, Max)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(4);
	node_int* b = new node_int(6);
	node_int* c = new node_int(5);
	t.appendNode(a);
	t.appendNode(b);
	t.appendNode(c);
	ASSERT_EQ(b, a->max());
}

TEST(Neighbours, Min)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(6);
	node_int* b = new node_int(2);
	node_int* c = new node_int(4);
	t.appendNode(a);
	t.appendNode(b);
	t.appendNode(c);
	ASSERT_EQ(b, a->min());
}

TEST(Neighbours, NextUp)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(6);
	node_int* b = new node_int(2);
	node_int* c = new node_int(4);
	t.insertAndRebalance(a);
	t.insertAndRebalance(b);
	t.insertAndRebalance(c);
	ASSERT_EQ(a, c->next());
	ASSERT_EQ(&t.m_end, a->next());
}

TEST(Neighbours, NextDown)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(2);
	node_int* b = new node_int(6);
	node_int* c = new node_int(4);
	t.insertAndRebalance(a);
	t.insertAndRebalance(b);
	t.insertAndRebalance(c);
	ASSERT_EQ(c, a->next());
	ASSERT_EQ(&t.m_end, b->next());
}

TEST(Neighbours, PrevUp)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(2);
	node_int* b = new node_int(6);
	node_int* c = new node_int(4);
	t.appendNode(a);
	t.appendNode(b);
	t.appendNode(c);
	ASSERT_EQ(a, c->previous());
	ASSERT_EQ(&t.m_end, a->previous());
}

TEST(Neighbours, PrevDown)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(6);
	node_int* b = new node_int(2);
	node_int* c = new node_int(4);
	t.appendNode(a);
	t.appendNode(b);
	t.appendNode(c);
	ASSERT_EQ(c, a->previous());
	ASSERT_EQ(&t.m_end, b->previous());
}

TEST(Deletion, not_invalidate_iterators)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	order_statistics_tree<int>::iterator its[10];
	for(int i = 0 ; i < 10; ++i)
	{
		its[i] = t.insert(tab[i]).first;
	}

	for(int j = 0 ; j < 10 ; ++j)
	{
		for(int i = j ; i < 10; ++i)
		{
			ASSERT_EQ(tab[i], *its[i]);
		}

		ASSERT_EQ(static_cast<size_t>(1), t.erase(*its[j]));
		ASSERT_EQ(static_cast<size_t>(9-j), t.size());
	}
}

TEST(Deletion, OneElement146)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}
	nodes[3] -> m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.deleteNode(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[3]);
	ASSERT_TRUE(nodes[1]->m_left == 0);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[3]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(4), t.size());
}

TEST(Deletion, OneElement2478)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}
	nodes[3] -> m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.deleteNode(nodes[3]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[4]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[1]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[4]);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(4), t.size());
}

TEST(Deletion, OneElement)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(6);
	t.appendNode(a);
	t.deleteNode(a);
	ASSERT_EQ(0, t.m_root);
	ASSERT_TRUE(t.empty());
	ASSERT_EQ(static_cast<size_t>(0), t.size());
}

TEST(Deletion, TwoElements12)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(6);
	node_int* b = new node_int(11);
	t.appendNode(a);
	t.appendNode(b);
	t.deleteNode(a);
	ASSERT_TRUE(t.m_root == b);
	ASSERT_EQ(0, b->m_left);
	ASSERT_EQ(0, b->m_right);
	ASSERT_TRUE(b->m_parent == &t.m_end);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(1), t.size());
}

TEST(Deletion, TwoElements21)
{
	order_statistics_tree<int> t;
	node_int* a = new node_int(11);
	node_int* b = new node_int(6);
	t.appendNode(a);
	t.appendNode(b);
	t.deleteNode(a);
	ASSERT_TRUE(t.m_root == b);
	ASSERT_EQ(0, b->m_left);
	ASSERT_EQ(0, b->m_right);
	ASSERT_TRUE(b->m_parent == &t.m_end);

	ASSERT_FALSE(t.empty());
	ASSERT_EQ(static_cast<size_t>(1), t.size());
}

TEST(DeletionFix, Case12)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}
	nodes[3] -> m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.RBdeleteFix(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[3]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false); // strange
	ASSERT_TRUE(nodes[2]->m_isRed == true); //strange
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);
}

TEST(DeletionFix, Case4)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}
	nodes[1] -> m_isRed = true;
	nodes[2] -> m_isRed = true;
	nodes[4] -> m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);

	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);

	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == true);

	t.RBdeleteFix(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[3]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);
}

TEST(DeletionFix, Case34)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}
	nodes[1] -> m_isRed = true;
	nodes[2] -> m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == true);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.RBdeleteFix(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[2]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == nodes[1]);
	ASSERT_TRUE(nodes[3]->m_left == 0);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == 0);
	ASSERT_TRUE(nodes[2]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[2]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);
}

TEST(DeletionFix, Case2a)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.RBdeleteFix(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[4]->m_isRed == false);
}

TEST(DeletionFix, Case2b)
{
	order_statistics_tree<int> t;
	int values[] = {0, 1, 2, 3, 4};
	node_int* nodes[5];
	for(int i = 0 ; i < 5 ; ++i)
	{
		nodes[i] = new node_int(values[i]);
	}

	nodes[1]->m_isRed = true;

	t.appendNode(nodes[1]);
	t.appendNode(nodes[3]);
	t.appendNode(nodes[0]);
	t.appendNode(nodes[4]);
	t.appendNode(nodes[2]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == true);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == false);
	ASSERT_TRUE(nodes[4]->m_isRed == false);

	t.RBdeleteFix(nodes[0]);

	ASSERT_TRUE(t.m_root == nodes[1]);
	ASSERT_TRUE(nodes[0]->m_left == 0);
	ASSERT_TRUE(nodes[1]->m_left == nodes[0]);
	ASSERT_TRUE(nodes[2]->m_left == 0);
	ASSERT_TRUE(nodes[3]->m_left == nodes[2]);
	ASSERT_TRUE(nodes[4]->m_left == 0);
	ASSERT_TRUE(nodes[0]->m_right == 0);
	ASSERT_TRUE(nodes[1]->m_right == nodes[3]);
	ASSERT_TRUE(nodes[2]->m_right == 0);
	ASSERT_TRUE(nodes[3]->m_right == nodes[4]);
	ASSERT_TRUE(nodes[4]->m_right == 0);
	ASSERT_TRUE(nodes[0]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[1]->m_parent == &t.m_end);
	ASSERT_TRUE(nodes[2]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[3]->m_parent == nodes[1]);
	ASSERT_TRUE(nodes[4]->m_parent == nodes[3]);
	ASSERT_TRUE(nodes[0]->m_isRed == false);
	ASSERT_TRUE(nodes[1]->m_isRed == false);
	ASSERT_TRUE(nodes[2]->m_isRed == false);
	ASSERT_TRUE(nodes[3]->m_isRed == true);
	ASSERT_TRUE(nodes[4]->m_isRed == false);
}

TEST(Iterators, begin_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::iterator it = t.begin();
	++it;
	ASSERT_EQ(2, *it); //get second element

	--it;
	ASSERT_EQ(1, *it); //back to begin

	ASSERT_THROW(--it, std::runtime_error ); // go before begin
}

TEST(Iterators, end_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::iterator it = t.end();
	--it;
	ASSERT_EQ(10, *it); //get last element

	++it;
	ASSERT_EQ(t.end(), it); //back to end

	ASSERT_THROW(++it, std::runtime_error ); // go after end
}

TEST(Iterators, incrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::iterator it = t.begin();
	ASSERT_EQ(1, *(it++));
	ASSERT_EQ(2, *it);

	it = t.begin();
	ASSERT_EQ(2, *(++it));
}

TEST(Iterators, decrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::iterator it = t.begin();
	++it;
	ASSERT_EQ(2, *(it--));
	ASSERT_EQ(1, *it);

	it = t.begin();
	++it;
	ASSERT_EQ(1, *(--it));
}

TEST(Iterators, position)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};

	for(int i = 0 ; i < 10; ++i)
	{
		t.insert(tab[i]);
	}

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t.begin(); it != t.end() ; ++it)
	{
		ASSERT_EQ(i, it.position());
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Iterators, loop)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t.begin(); it != t.end() ; ++it)
	{
		ASSERT_EQ(*it, tab[i]);
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Citerators, begin_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_iterator it = t.cbegin();
	++it;
	ASSERT_EQ(2, *it); //get second element

	--it;
	ASSERT_EQ(1, *it); //back to begin

	ASSERT_THROW(--it, std::runtime_error ); // go before begin
}

TEST(Citerators, end_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_iterator it = t.cend();
	--it;
	ASSERT_EQ(10, *it); //get last element

	++it;
	ASSERT_EQ(t.cend(), it); //back to end

	ASSERT_THROW(++it, std::runtime_error ); // go after end
}

TEST(Citerators, incrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_iterator it = t.cbegin();
	ASSERT_EQ(1, *(it++));
	ASSERT_EQ(2, *it);

	it = t.cbegin();
	ASSERT_EQ(2, *(++it));
}

TEST(Citerators, decrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_iterator it = t.cbegin();
	++it;
	ASSERT_EQ(2, *(it--));
	ASSERT_EQ(1, *it);

	it = t.cbegin();
	++it;
	ASSERT_EQ(1, *(--it));
}

TEST(Citerators, position)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};

	for(int i = 0 ; i < 10; ++i)
	{
		t.insert(tab[i]);
	}

	size_t i = 0;
	for(order_statistics_tree<int>::const_iterator it = t.cbegin(); it != t.cend() ; ++it)
	{
		ASSERT_EQ(i, it.position());
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Citerators, loop)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);

	size_t i = 0;
	for(order_statistics_tree<int>::const_iterator it = t.cbegin(); it != t.cend() ; ++it)
	{
		ASSERT_EQ(*it, tab[i]);
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Riterators, begin_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::reverse_iterator it = t.rbegin();
	++it;
	ASSERT_EQ(9, *it); //get second element

	--it;
	ASSERT_EQ(10, *it); //back to begin

	ASSERT_THROW(--it, std::runtime_error ); // go before begin
}

TEST(Riterators, end_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::reverse_iterator it = t.rend();
	--it;
	ASSERT_EQ(1, *it); //get last element

	++it;
	ASSERT_EQ(t.rend(), it); //back to end

	ASSERT_THROW(++it, std::runtime_error ); // go after end
}

TEST(Riterators, incrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(size_t i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::reverse_iterator it = t.rbegin();
	ASSERT_EQ(10, *(it++));
	ASSERT_EQ(9, *it);

	it = t.rbegin();
	ASSERT_EQ(9, *(++it));
}

TEST(Riterators, decrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::reverse_iterator it = t.rbegin();
	++it;
	ASSERT_EQ(9, *(it--));
	ASSERT_EQ(10, *it);

	it = t.rbegin();
	++it;
	ASSERT_EQ(10, *(--it));
}

TEST(Riterators, position)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};

	for(int i = 0 ; i < 10; ++i)
	{
		t.insert(tab[i]);
	}

	size_t i = 0;
	for(order_statistics_tree<int>::reverse_iterator it = t.rbegin(); it != t.rend() ; ++it)
	{
		ASSERT_EQ(i, it.position());
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Riterators, loop)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);
	std::reverse(tab, tab+10);

	size_t i = 0;
	for(order_statistics_tree<int>::reverse_iterator it = t.rbegin(); it != t.rend() ; ++it)
	{
		ASSERT_EQ(*it, tab[i]);
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(CRiterators, begin_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_reverse_iterator it = t.crbegin();
	++it;
	ASSERT_EQ(9, *it); //get second element

	--it;
	ASSERT_EQ(10, *it); //back to begin

	ASSERT_THROW(--it, std::runtime_error ); // go before begin
}

TEST(CRiterators, end_iterator_arithmetic)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_reverse_iterator it = t.crend();
	--it;
	ASSERT_EQ(1, *it); //get last element

	++it;
	ASSERT_EQ(t.crend(), it); //back to end

	ASSERT_THROW(++it, std::runtime_error ); // go after end
}

TEST(CRiterators, incrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_reverse_iterator it = t.crbegin();
	ASSERT_EQ(10, *(it++));
	ASSERT_EQ(9, *it);

	it = t.crbegin();
	ASSERT_EQ(9, *(++it));
}

TEST(CRiterators, decrementation)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(int i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	order_statistics_tree<int>::const_reverse_iterator it = t.crbegin();
	++it;
	ASSERT_EQ(9, *(it--));
	ASSERT_EQ(10, *it);

	it = t.crbegin();
	++it;
	ASSERT_EQ(10, *(--it));
}

TEST(CRiterators, position)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};

	for(int i = 0 ; i < 10; ++i)
	{
		t.insert(tab[i]);
	}

	size_t i = 0;
	for(order_statistics_tree<int>::const_reverse_iterator it = t.crbegin(); it != t.crend() ; ++it)
	{
		ASSERT_EQ(i, it.position());
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(CRiterators, loop)
{
	order_statistics_tree<int> t;
	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(size_t i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);
	std::reverse(tab, tab+10);

	size_t i = 0;
	for(order_statistics_tree<int>::const_reverse_iterator it = t.crbegin(); it != t.crend() ; ++it)
	{
		ASSERT_EQ(*it, tab[i]);
		++i;
	}

	ASSERT_EQ(t.size(), i);
}

TEST(Operators, position)
{
	order_statistics_tree<int> t;
	ASSERT_EQ(t.end(), t[0]);

	int tab[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	node_int* ptrs[10];
	for(size_t i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);

	for(size_t i = 0 ; i < 10; ++i)
	{
		ASSERT_EQ(tab[i], (*t[i]) );
	}

	ASSERT_EQ(t.end(), t[10000000]);
}

TEST(Operators, position_randomInsertAndErase)
{
	srand ( static_cast<unsigned int>(time(0)) );
	order_statistics_tree<int> t;
	std::vector<int> v;

	const size_t SIZE = 5000;
	const size_t REPLACE_SIZE = 100;

	for(size_t i = 0 ; i < SIZE ; ++i)
	{
		int tmp = rand() % 1000000;
		t.insert(tmp);
		v.push_back(tmp);
	}

	std::sort(v.begin(), v.end());

	for(int tests = 0 ; tests < 2000 ; ++tests)
	{
		for(size_t liczby = 0 ; liczby < REPLACE_SIZE ; ++liczby)
		{
			int indeks = rand() % (SIZE-liczby);
			t.erase(t[indeks]);
			v.erase(v.begin() + indeks);
		}

		for(size_t liczby = 0 ; liczby < REPLACE_SIZE ; ++liczby)
		{
			int val = rand() % 1000000;
			t.insert(val);
			v.push_back(val);
		}

		std::sort(v.begin(), v.end());

		for(size_t j = 0 ; j < SIZE ; ++j)
		{
			ASSERT_EQ(v[j], *t[j]);
		}
	}
}

node_int* ptrs[10];
TEST(Operators, position_Rem_Add)
{
	order_statistics_tree<int> t;
	int tab[10] = {61, 91, 21, 11, 41, 51, 71, 101, 81, 31};
	
	for(size_t i = 0 ; i < 10; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab, tab+10);
	std::sort(ptrs, ptrs+10, NodeCmp() );

	for(size_t i = 0 ; i < 10; ++i)
	{
		ASSERT_EQ(tab[i], *t[i] );
	}

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		for(size_t j = 0 ; j < 10 - i; ++j)
		{
			ASSERT_EQ(tab[i+j], *t[j] );
		}
		t.deleteNode(ptrs[i]);
	}
	ASSERT_EQ(0, t.m_root);
}

TEST(Operators, position_Random)
{
	order_statistics_tree<int> t;
	const size_t SIZE = 10;
	srand ( static_cast<unsigned int>(time(0)) );
	for(int testy = 0 ; testy < 100000 ; ++testy)
	{
	std::vector<int> tab(SIZE);
	for(size_t i = 0 ; i < SIZE ; ++i)
	{
		tab[i] = rand()%100;
	}

	for(size_t i = 0 ; i < SIZE; ++i)
	{
		ptrs[i] = new node_int( tab[i] );
		t.insertAndRebalance(ptrs[i]);
	}

	std::sort(tab.begin(), tab.end());
	std::sort(ptrs, ptrs+SIZE, NodeCmp() );

	for(size_t i = 0 ; i < SIZE; ++i)
	{
		ASSERT_EQ(tab[i], *t[i] );
	}

	for(size_t i = 0 ; i < SIZE ; ++i)
	{
		for(size_t j = 0 ; j < SIZE - i; ++j)
		{
			ASSERT_EQ(tab[i+j], *t[j] );
		}
		t.deleteNode(ptrs[i]);
	}
	ASSERT_EQ(0, t.m_root);
	}
}

TEST(Miscellaneous, TreeForString_basic)
{
	order_statistics_tree<std::string, LengthyStringComparator> t;

	std::string tab[10] = { "AAA", "AAB", "BA", "BBB", "C", "D", "CDD", "DA", "A", ""};
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		t.insert(tab[i]);
	}

	ASSERT_EQ(static_cast<size_t>(10), t.size());

	std::sort(tab, tab+10, LengthyStringComparator() );
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		ASSERT_EQ(tab[i], *t[i]);
	}
}

TEST(Miscellaneous, TreeForString_add_and_delete)
{
	order_statistics_tree<std::string, LengthyStringComparator> t;

	std::string tab[10] = { "AAA", "AAB", "BA", "BBB", "C", "D", "CDD", "DA", "A", ""};
	order_statistics_tree<std::string, LengthyStringComparator>::iterator its[10];
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		its[i] = t.insert(tab[i]).first;
	}

	std::sort(tab, tab+10, LengthyStringComparator());
	std::sort(its, its+10, IteratorCmp<std::string, LengthyStringComparator>());

	ASSERT_EQ(static_cast<size_t>(10), t.size());

	t.erase(tab[0]);
	ASSERT_EQ(static_cast<size_t>(9), t.size());

	ASSERT_EQ(*its[2], *(++its[1]));
	ASSERT_EQ(its[6], t.find(tab[6]));

	ASSERT_EQ(its[5], t.equal_range(*its[5]).first);
	ASSERT_EQ(its[6], t.equal_range(*its[5]).second);

}

bool redness;

template <typename ValueType>
std::pair<int, int> validator(typename order_statistics_tree<ValueType>::node_pointer l_root, int height, bool red)
{
	if(l_root == 0)
	{
		return std::make_pair(height, height);
	}

	redness |= l_root->m_isRed && red;

	std::pair<int, int> left = validator<ValueType>(l_root->m_left, height + 1, l_root->m_isRed);
	std::pair<int, int> right = validator<ValueType>(l_root->m_right, height + 1, l_root->m_isRed);

	int tab[4];
	tab[0] = left.first;
	tab[1] = left.second;
	tab[2] = right.first;
	tab[3] = right.second;

	std::sort(tab, tab+4);

	return std::make_pair(tab[0], tab[3]);
}

TEST(Miscellaneous, TreeStructureValidator_basic)
{
	for(size_t j = 0 ; j < 100 ; ++j)
	{
		redness = false;
		order_statistics_tree<int> t;
		srand ( static_cast<unsigned int>(time(0)) );
		for(size_t i = 0 ; i < 100000 ; ++i)
		{
			t.insert(rand() % 100000000 );
		}

		std::pair<int, int> res = validator<int>(t.m_root, 0, true);
		ASSERT_TRUE( 2 * res.first >= res.second);
		ASSERT_FALSE(redness);
	}
}

TEST(Miscellaneous, CompatibilityWithSTL)
{
	order_statistics_tree<std::string, LengthyStringComparator> t;

	std::string tab[10] = { "AAA", "AAB", "BA", "BBB", "C", "D", "CDD", "DA", "A", ""};
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		t.insert(tab[i]);
	}

	ASSERT_EQ(static_cast<size_t>(10), t.size());

	std::sort(tab, tab+10, LengthyStringComparator() );
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		ASSERT_EQ(tab[i], *t[i]);
	}

	std::set<std::string> s1(t.begin(), t.end());
	ASSERT_EQ(static_cast<size_t>(10), s1.size());

	std::set<std::string> s2(t.cbegin(), t.cend());
	ASSERT_EQ(static_cast<size_t>(10), s2.size());

	std::set<std::string> s3(t.rbegin(), t.rend());
	ASSERT_EQ(static_cast<size_t>(10), s3.size());

	std::set<std::string> s4(t.crbegin(), t.crend());
	ASSERT_EQ(static_cast<size_t>(10), s4.size());

	order_statistics_tree<std::string, LengthyStringComparator>::iterator it = t.begin();
	order_statistics_tree<std::string, LengthyStringComparator>::iterator it2 = t.end();
	std::swap(it, it2);
	ASSERT_EQ(t.begin(), it2);
	ASSERT_EQ(t.end(), it);
}

struct Pair
{
	int key;
	int extraVal;

	Pair(int key = 0, int extraVal = 0) : key(key), extraVal(extraVal) {}

	bool operator<(const Pair& rhs) const
	{
		return key < rhs.key;
	}
};

TEST(Miscellaneous, Multiinsertion_Stability)
{
	order_statistics_tree<Pair> t;

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		t.insert(Pair(1, i));
	}

	for(int i = 0 ; i < 10 ; ++i)
	{
		ASSERT_EQ(i, t[i]->extraVal);
	}
}

TEST(Miscellaneous, EmptyTree)
{
	order_statistics_tree<std::string> t;
	ASSERT_EQ(t.end(), t.find("A"));
}

TEST(Bounds, LowerBound)
{
	order_statistics_tree<int> t;
	int values[] = {15, 15, 15, 15, 18};
	order_statistics_tree<int>::iterator its[5];
	for(size_t i = 0 ; i < 5 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(its[0], t.lower_bound(14));
	ASSERT_EQ(its[0], t.lower_bound(15));
	ASSERT_EQ(its[4], t.lower_bound(16));
	ASSERT_EQ(t.end(), t.lower_bound(19));
}

TEST(Bounds, UpperBound)
{
	order_statistics_tree<int> t;
	int values[] = {15, 15, 15, 15, 18};
	order_statistics_tree<int>::iterator its[5];
	for(size_t i = 0 ; i < 5 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(its[0], t.upper_bound(14));
	ASSERT_EQ(its[4], t.upper_bound(15));
	ASSERT_EQ(its[4], t.upper_bound(16));
	ASSERT_EQ(t.end(), t.upper_bound(19));
}

TEST(Bounds, EqualRange)
{
	order_statistics_tree<int> t;
	int values[] = {15, 15, 15, 15, 18};
	order_statistics_tree<int>::iterator its[5];
	for(size_t i = 0 ; i < 5 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(its[0], t.equal_range(14).first);
	ASSERT_EQ(its[0], t.equal_range(14).second);
	ASSERT_EQ(its[0], t.equal_range(15).first);
	ASSERT_EQ(its[4], t.equal_range(15).second);
	ASSERT_EQ(its[4], t.equal_range(16).first);
	ASSERT_EQ(its[4], t.equal_range(16).second);
	ASSERT_EQ(t.end(), t.equal_range(19).first);
	ASSERT_EQ(t.end(), t.equal_range(19).second);
}

TEST(Bounds, Find)
{
	order_statistics_tree<int> t;
	int values[] = {15, 15, 15, 15, 18};
	order_statistics_tree<int>::iterator its[5];
	for(size_t i = 0 ; i < 5 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(t.end(), t.find(14));
	ASSERT_EQ(its[0], t.find(15));
	ASSERT_EQ(t.end(), t.find(16));
	ASSERT_EQ(t.end(), t.find(19));
}

TEST(Bounds, Count)
{
	order_statistics_tree<int> t;
	int values[] = {15, 15, 15, 15, 18};
	order_statistics_tree<int>::iterator its[5];
	for(size_t i = 0 ; i < 5 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(static_cast<size_t>(0), t.count(14));
	ASSERT_EQ(static_cast<size_t>(4), t.count(15));
	ASSERT_EQ(static_cast<size_t>(0), t.count(16));
	ASSERT_EQ(static_cast<size_t>(0), t.count(19));
	ASSERT_EQ(static_cast<size_t>(1), t.count(18));
}

TEST(Functions, Clear)
{
	order_statistics_tree<int> t;
	int values[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	order_statistics_tree<int>::iterator its[10];
	for(size_t i = 0 ; i < 10 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(static_cast<size_t>(10), t.size());
	t.clear();
	ASSERT_EQ(static_cast<size_t>(0), t.size());
	ASSERT_EQ(0, t.m_root);
	ASSERT_EQ(0, t.m_end.m_parent);
	ASSERT_EQ(0, t.m_end.m_right);
	ASSERT_EQ(0, t.m_end.m_left);
}

TEST(Functions, Swap)
{
	order_statistics_tree<int> t[2];
	int values[2][10] = {{6, 9, 2, 1, 4, 5, 7, 10, 8, 3}, {600, 900, 200, 100, 400, 500, 700, 1000, 800, 300}};
	order_statistics_tree<int>::iterator its[2][10];

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		its[0][i] = t[0].insert(values[0][i]).first;
		its[1][i] = t[1].insert(values[1][i]).first;
	}

	ASSERT_EQ(static_cast<size_t>(10), t[0].size());
	ASSERT_EQ(static_cast<size_t>(10), t[1].size());

	std::sort(its[0], its[0]+10, IteratorCmp<int>());
	std::sort(its[1], its[1]+10, IteratorCmp<int>());

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t[0].begin() ; it != t[0].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[0][i], it);
	}

	i = 0;
	for(order_statistics_tree<int>::iterator it = t[1].begin() ; it != t[1].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[1][i], it);
	}

	t[0].swap(t[1]);

	i = 0;
	for(order_statistics_tree<int>::iterator it = t[0].begin() ; it != t[0].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[1][i], it);
	}

	i = 0;
	for(order_statistics_tree<int>::iterator it = t[1].begin() ; it != t[1].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[0][i], it);
	}
}

TEST(Functions, CopyConstructor)
{
	order_statistics_tree<int> t;
	int values[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	order_statistics_tree<int>::iterator its[10];

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		its[i] = t.insert(values[i]).first;
	}

	ASSERT_EQ(static_cast<size_t>(10), t.size());

	std::sort(its, its+10, IteratorCmp<int>());

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t.begin() ; it != t.end() ; ++it, ++i)
	{
		ASSERT_EQ(its[i], it);
	}

	order_statistics_tree<int> copy(t);

	i = 0;
	for(order_statistics_tree<int>::iterator it = t.begin() ; it != t.end() ; ++it, ++i)
	{
		ASSERT_EQ(its[i], it);
	}

	i = 0;
	for(order_statistics_tree<int>::iterator it = copy.begin() ; it != copy.end() ; ++it, ++i)
	{
		ASSERT_EQ(*its[i], *it);
	}
}

TEST(Functions, ConstructorFromIterators)
{
	std::set<int> s;
	int values[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		s.insert(values[i]);
	}

	ASSERT_EQ(10U, s.size());
	std::sort(values, values+10);

	order_statistics_tree<int> t(s.begin(), s.end());
	ASSERT_EQ(10U, t.size());

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t.begin() ; it != t.end() ; ++it, ++i)
	{
		ASSERT_EQ(values[i], *it);
	}
}

TEST(Functions, AssignOperator)
{
	order_statistics_tree<int> t[2];
	int values[10] = {6, 9, 2, 1, 4, 5, 7, 10, 8, 3};
	order_statistics_tree<int>::iterator its[10];

	for(size_t i = 0 ; i < 10 ; ++i)
	{
		its[i] = t[0].insert(values[i]).first;
	}

	ASSERT_EQ(10U, t[0].size());

	std::sort(its, its+10, IteratorCmp<int>());

	size_t i = 0;
	for(order_statistics_tree<int>::iterator it = t[0].begin() ; it != t[0].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[i], it);
	}

	t[1] = t[0];

	i = 0;
	for(order_statistics_tree<int>::iterator it = t[0].begin() ; it != t[0].end() ; ++it, ++i)
	{
		ASSERT_EQ(its[i], it);
	}

	i = 0;
	for(order_statistics_tree<int>::iterator it = t[1].begin() ; it != t[1].end() ; ++it, ++i)
	{
		ASSERT_EQ(*its[i], *it);
	}
}

TEST(Functions, Erase)
{
	order_statistics_tree<int> t;
	for(size_t i = 0 ; i < 20 ; ++i)
	{
		t.insert(i);
		t.insert(100);
	}

	ASSERT_EQ(40U, t.size());
	ASSERT_EQ(20U, t.erase(100));

	t.erase(t[0], t[15]);
	ASSERT_EQ(5U, t.size());

	t.erase(18);
	ASSERT_EQ(4U, t.size());

	t.erase(50);
	ASSERT_EQ(4U, t.size());

	t.erase(t[100000], t[20000]);
	ASSERT_EQ(4U, t.size());
}

template <
	class PropObject,
	class PropType, 
	PropType PropObject::* Prop
	>
struct FieldRef
{
	 typedef PropObject argument_type;
	 typedef PropType result_type;

	 typedef PropType value_type;

	 PropType operator()( const PropObject& p ) const { return p.*Prop; }
	 PropType& operator()( PropObject& p ) const { return p.*Prop; }
};

template <typename T1, typename T2>
std::ostream& operator<<( std::ostream& s, const std::pair< T1, T2>& p ) {
	s << "(" << p.first << "," << p.second << ")";
	return s;
}


TEST(RangeUpdate, Add)
{
	typedef double key_t;
	typedef int weight_t;
	typedef std::pair< key_t, weight_t > pair_t;
	typedef FieldRef< pair_t, key_t, &pair_t::first > GetFirst;
	typedef FieldRef< pair_t, weight_t, &pair_t::second > GetSecond;

	typedef order_statistics_tree< pair_t, std::less<key_t>, GetFirst, GetSecond > ost_t;
	
	ost_t t;
	for(size_t i = 0 ; i < 20 ; ++i)
		 t.insert( std::make_pair( static_cast<key_t>(i), static_cast<weight_t>( 2 * i) ) );

	t.insert( std::make_pair( static_cast<key_t>(5.5), static_cast<weight_t>( 239) ) );
	
	t.addWeight( t.find( 1 ), t.find( 5 ), 100 );

	ost_t::iterator ub5 = t.upper_bound( static_cast<key_t>( 5 ) );
	ost_t::iterator ub4 = t.upper_bound( static_cast<key_t>( 4 ) );
	
	// std::cerr << "upper bound: " << *ub5 << std::endl;
	// std::cerr << "upper bound: " << *ub4 << std::endl;
	// std::cerr << "less? " << t.is_less( ub4, ub5 ) << std::endl;
	
	// for(ost_t::const_iterator i = t.begin(); i != t.end(); i++ )
	// 	 std::cerr << "at " << i.position() << " value=" << *i << std::endl;

	ASSERT_EQ( t.find( static_cast<key_t>( 1 ) )->second, static_cast<weight_t>( 102 ) );
	
}

TEST(EqualVals, EqVal1)
{
	typedef order_statistics_tree< double > ost_t;
	typedef ost_t::iterator ost_iter_t;
	
	ost_t t;

	ost_iter_t i0 = t.insert( 2.0 ).first;
	ost_iter_t i1 = t.insert( 2.0 ).first;
	ost_iter_t i2 = t.insert( 2.0 ).first;

	std::cerr << "i0=" << i0.position() << " i1=" << i1.position() << " i2=" << i2.position() << std::endl;
	
}

TEST(PSumSelect, PSumSelect1)
{
	using std::pair;
	using std::make_pair;
	using std::cerr;
	using std::endl;

	typedef double key_t;
	typedef int weight_t;
	typedef pair< key_t, weight_t > pair_t;
	typedef FieldRef< pair_t, key_t, &pair_t::first > GetFirst;
	typedef FieldRef< pair_t, weight_t, &pair_t::second > GetSecond;

	typedef order_statistics_tree< pair_t, std::less<key_t>, GetFirst, GetSecond > ost_t;

	ost_t t;
	t.insert( make_pair( .2, 5 ) );
	t.insert( make_pair( .3, 7 ) );
	t.insert( make_pair( .5, 3 ) );
	for ( int i = 0; i <= 15; i++ ) {
		int residue = -10000;
		ost_t::const_iterator it = t.upper_bound_for_inclusive_partial_sum( i, &residue );
		cerr << "i=" << i << " it.pos=" << it.position() << " residue=" << residue;
		if ( it == ost_t::const_iterator( t.end() ) ) cerr << " it=END";
		else cerr << " first=" << it->first << " second=" << it->second;
		cerr << endl;
	}
}

int main(int argc, char **argv)
{
//	::testing::GTEST_FLAG(filter) = "*add_and_delete*";
  ::testing::InitGoogleTest(&argc, argv);

  int testResult = RUN_ALL_TESTS();

	std::cerr << "testResult=" << testResult << std::endl;

  return testResult;
}
