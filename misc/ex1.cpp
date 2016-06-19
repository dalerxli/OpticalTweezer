#include <iostream>

using namespace std;

struct Node
{
				int data;
				Node *next;
};


//1st node:
void firstNode(struct Node *head,int n)
{
				head -> data = n;
				head -> next = NULL;
				cout << "Adresse: " << &head << "\n";
				cout << "Wert: " << head << "\n";
}

//appending
void addNode(struct Node *head,int n)
{
				Node *nextNode = new Node;
				nextNode -> data = n;
				nextNode -> next = NULL;

				Node *cur = head;
				while(cur)
				{
								if(cur -> next == NULL)
								{
												cur -> next = nextNode;
												return;
								}
								cur = cur -> next;
				}
}

void insertFront(struct Node **head,int n)
{
				struct Node *newHead = new Node;
				newHead -> data = n;
				newHead -> next = *head;
				*head = newHead;
}

struct Node *insFront(struct Node **head,int n)
{
				struct Node *temp = new Node;
				temp -> data = n;
				temp -> next = *head;
				*head = temp;
				return *head;
}

struct Node *searchNode(struct Node *head,int n)
{
		Node *current = head;
		while(current)
		{
						if(current -> data == n)
										return current;
						current = current -> next;
		}
		cout << "No Node " << n << " in list.\n";
}

bool deleteNode(struct Node **head, Node *toDelete)
{
				Node *current = *head;
				if(toDelete == *head)
				{
								*head = current -> next;
								delete toDelete;
								return true;
				}

				while(current)
				{
								if(current -> next == toDelete)
								{
												current -> next = toDelete -> next;
												delete toDelete;
												return true;
								}
								current = current -> next;
				}
				return false;
}


//reverse the list
struct Node *reverseList(struct Node **head)
{
				Node *parent = *head;
				Node *me = parent -> next;
				Node *child = me -> next;

				parent -> next = NULL; // parent is last item 
				while(child)
				{
								me -> next = parent;
								parent = me;
								me = child;
								child = child -> next;
				}
				me -> next = parent;
				*head = me;
				return *head;
}

//copy a list
void copyLinkedList(struct Node *node, struct Node **newList)
{
				if(node != NULL)
				{
								*newList = new Node;
								(*newList) -> data = node -> data;
								(*newList) -> next = NULL;
								copyLinkedList(node->next,&((*newList)->next));
				}
}

//compare lists
//0 ... lists are different
//1 ... lists are the same
int compareLists(struct Node *node1, struct Node *node2)
{
				static int flag;

				if(node1 == NULL && node2 == NULL)
								flag = 1;
				else
				{
								if(node1 == NULL || node2 == NULL)
												flag = 0;
								else if(node1->data != node2 -> data)
												flag = 0;
								else
												compareLists(node1->next,node2->next);
				}

				return flag;
				
}


void deleteList(struct Node **node)
{
				struct Node *temporary;
				while(*node)
				{
								temporary = *node;
								*node = temporary -> next;
								delete temporary;
				}
}


void displayList(struct Node *head)
{
				Node *list = head;
				while(list)
				{
								cout << list -> data << " ";
								list = list -> next;
				}

				cout << "\n\n";
}


int main()
{
				struct Node *newHead;
				struct Node *head = new Node;

				cout << "Adresse: " << &head << "\n";
				cout << "Wert: " << head << "\n";

				firstNode(head,10);
				displayList(head);

				addNode(head,20);
				displayList(head);


				addNode(head,30);
				displayList(head);


				addNode(head,35);
				displayList(head);


				addNode(head,40);
				displayList(head);


				insertFront(&head,5);
				displayList(head);

				//head = insFront(&head,5);
				//displayList(head);
				
				//addNode(head,5);
				//displayList(head);	
				
				struct Node *iter = head;
				int numDel = 5;

				while(iter)
				{
								if(iter -> data == numDel)
								{
												deleteNode(&head,iter);
								}
								iter = iter -> next;
				}
				displayList(head);
				/*
				 *Node *ptrDelete = searchNode(head,numDel);
				 *if(deleteNode(&head,ptrDelete))
				 *        cout << "Node " << numDel << " deleted!\n";
				 *displayList(head);
				 */

				/*
				 *cout << "Reversed List: \n";
				 *reverseList(&head);
				 *displayList(head);
				 */

				/*
				 *cout << "Copying list...\n";
				 *copyLinkedList(head,&newHead);
				 *displayList(newHead);
				 */

				return 0;
}





























