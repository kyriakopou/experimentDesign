/*
 *  dynspace.h
 *  Attraction
 *
 *  Created by David Spieler on 2/17/12.
 *  Copyright 2012 Saarland University. All rights reserved.
 *
 */

#ifndef __DYNSPACE_H__
#define __DYNSPACE_H__

#include "common.h"

namespace attraction
{
	template <typename K, typename V, typename Equals, typename Hash, uint MaxFanout>
	class DynSpace
	{
	public:
		struct Entry
		{
			K	key;		/* key */
			V	value;		/* value */
			
			uint	fanout;        	/* #successor states = #possible reactions, INVALID if still unexplored */
			real	tp[MaxFanout]; 	/* transition rates to the successor state*/
			int		tn[MaxFanout];	/*transition number*/
			real	tc[MaxFanout];	/*transition constant*/
			Entry*	s[MaxFanout];  	/* successor entries */
			real	exitRate;		/* exit rate */
			real	iexitRate;  	/* 1/exitRate */
			
			Entry*	h_next;		/* next pointer in bucket */
			
			Entry()
			{	
				
				for (uint i=0; i < MaxFanout; i++)
				{
					tp[i] = 0;
					s[i] = 0;
				}
				fanout = INVALID;
				exitRate = 0;
				iexitRate = 0;
				h_next = NULL;


			}
		};
	protected:
		Entry*		m_hash_table;
		Entry*		m_storage;
		
		unsigned int	m_hash_table_size;
		unsigned int	m_storage_size;
		unsigned int	m_nr_states;
		unsigned int	m_nr_hash_entries;
		unsigned int	m_next_free;
		unsigned int	m_pos_table;
		Entry*		m_unused;
		Entry*		m_pos_storage;
		
	public:
		DynSpace()
		{
			m_hash_table_size = 0;
			m_storage_size    = 0;
			m_hash_table      = 0;
			m_storage         = 0;
			m_next_free       = 0;
			m_unused          = 0;
			m_nr_states       = 0;
			m_nr_hash_entries = 0;
			m_pos_table       = 0;
			m_pos_storage     = 0;
		}
		
		~DynSpace()
		{
			clear();
		}
		
		void setup(unsigned int hash_table_size, unsigned int storage_size)
		{
			clear();
			
			m_hash_table_size = hash_table_size;
			m_storage_size    = storage_size;
			
			m_hash_table = new Entry[m_hash_table_size];
			m_storage    = new Entry[m_storage_size];
			m_unused     = &(m_storage[m_storage_size-1]);
			
			for (unsigned int i = 0; i < m_hash_table_size; i++)
			{
				m_hash_table[i].h_next = m_unused;
			}
			for (unsigned int i = 0; i < m_storage_size; i++)
			{
				m_storage[i].h_next = m_unused;
			}
		}
		
		void clear()
		{
			if (m_hash_table)
				delete [] m_hash_table;
			if (m_storage)
				delete [] m_storage;
			
			m_hash_table_size = 0;
			m_storage_size    = 0;
			m_hash_table      = 0;
			m_storage         = 0;
			m_next_free       = 0;
			m_unused          = 0;
			m_nr_states       = 0;
			m_nr_hash_entries = 0;
			m_pos_table       = 0;
			m_pos_storage     = 0;
		}
		
		/*hashmap: state->entry_ptr
		const& is used for speed*/
		Entry* lookup(const K& key)
		{
			unsigned int position = (Hash() (key)) % m_hash_table_size;
			
			if (m_hash_table[position].h_next == m_unused)
				return 0;
			
			Entry* entry = &(m_hash_table[position]);
			while (entry && !(Equals() (entry->key, key)))
				entry = entry->h_next;
			
			return entry;
		}
		
		Entry* add(const K& key)
		{
			unsigned position = (Hash() (key)) % m_hash_table_size;
			
			Entry* entry = &(m_hash_table[position]);
			
			if (entry->h_next == m_unused) 
			{
				entry->key      = key;
				entry->h_next   = 0;
				
				m_nr_states++;
				m_nr_hash_entries++;
				
				return entry;
			} else {
				while (entry->h_next)
					entry = entry->h_next;
				
				Entry *newEntry = &(m_storage[m_next_free++]);
				
				//assert(newEntry != m_unused);
				
				newEntry->key      = key;
				newEntry->h_next   = 0;
				
				entry->h_next = newEntry;
				
				m_nr_states++;
				
				return newEntry;
			}
		}
		
		Entry* ensureExistence(const K& key)
		{
			unsigned position = (Hash() (key)) % m_hash_table_size;
			
			Entry* entry = &(m_hash_table[position]);
			
			if (entry->h_next == m_unused) 
			{
				entry->key      = key;
				entry->h_next   = 0;
				
				m_nr_states++;
				m_nr_hash_entries++;
				
				return entry;
			} else {
				if (Equals() (entry->key, key))
					return entry;
				
				while (entry->h_next)
				{
					entry = entry->h_next;
					if (Equals() (entry->key, key))
						return entry;
				}
				
				Entry* newEntry = &(m_storage[m_next_free++]);
				
				//assert(newEntry != m_unused);
				
				newEntry->key      = key;
				newEntry->h_next   = 0;
				
				entry->h_next = newEntry;
				
				m_nr_states++;
				
				return newEntry;
			}
		}
	
		unsigned int get_hash_table_size() { return m_hash_table_size; };
		unsigned int get_storage_size()    { return m_storage_size; };
		unsigned int get_nr_States()       { return m_nr_states; };
		unsigned int get_nr_hash_entries() { return m_nr_hash_entries; };		
	};
}

#endif
