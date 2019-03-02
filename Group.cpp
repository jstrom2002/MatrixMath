#pragma once
#include "stdafx.h"

Group::Group() { GroupPrecision = 0; }

Group::~Group() {}

Group::Group(int n) {
		for (int i = 0; i < n; ++i) {
			element.push_back(i + 1);
		}
		GroupPrecision = 0;
	}

Group::Group(std::vector<int> v) {
		GroupPrecision = 0;
		if (v.size() < 1) { return; }
		if (v.size() == 1) {
			int n = v[0];
			for (int i = 0; i < n; ++i) {
				element.push_back(i + 1);
			}
			return;
		}
		if (v.size() > 1) {
			for (int i = 0; i < v.size(); ++i) {
				element.push_back(v[i]);
			}
		}
	}


	double Group::order() { return factorial((double)element.size()); }

	void Group::identity() {
		for (auto i = element.begin(); i != element.end(); ++i) {
			std::rotate(std::upper_bound(element.begin(), i, *i), i, i + 1);
		}
	}

	void Group::swap(int a, int b) {
		std::swap(element.begin() + a, element.begin() + b);
	}

	void Group::cycleLeft() {
		std::rotate(element.begin(), element.begin() + 1, element.end());
	}

	void Group::cycleRight() {
		std::rotate(element.rbegin(), element.rbegin() + 1, element.rend());
	}

	void Group::groupActionRight(std::vector<int> act) {
		//Cycle notation -- i.e (1,2,3)*(1,2) where (1,2) is the group action.		
		std::vector<int> pos;
		for (int i = 0; i < act.size(); ++i) {
			pos.push_back(findPosition(element, act[i]));
		}
		std::vector<int> copy = element;

		for (int i = 0; i < pos.size() - 1; ++i) {
			copy[pos[i + 1]] = element[pos[i]];
		}
		copy[pos[0]] = element[pos[pos.size() - 1]];
		element = copy;
	}

	void Group::groupActionLeft(std::vector<int> act) {
		//Cycle notation -- i.e. (1,2)*(1,2,3) = (1,2,3) where (1,2) is the group action
		std::vector<int> pos;
		for (int i = 0; i < act.size(); ++i) {
			pos.push_back(findPosition(element, act[i]));
		}
		std::vector<int> copy = element;

		for (int i = pos.size() - 2; i >= 0; --i) {
			copy[pos[i]] = element[pos[i + 1]];
		}
		copy[pos[pos.size() - 1]] = element[pos[0]];
		element = copy;
	}

	int Group::get(int i) { return element[i]; }
	void Group::set(int i, int x) { element[i] = x; }

	Vector Group::toVector() { return Vector(STLIntToDouble(element)); }

	std::wstring Group::toString() {
		std::wostringstream sstr;
		sstr << L"(";
		for (int i = 0; i < element.size(); ++i) {
			if (element[i]<999999999999999999 && element[i]> -999999999999999999) {
				sstr << std::fixed << std::setprecision(0) << element[i];
				if (i < element.size() - 1) { sstr << L","; }
			}
		}
		sstr << L")";
		std::wstring answer = sstr.str();
		sstr.clear();
		return answer;
	}

	void Group::display() { std::wcout << toString() << L"\n"; }